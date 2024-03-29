import itertools
import os
import sys
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
from functools import partial
from pathlib import Path
from timeit import default_timer as timer
from spapros.selection.selection_methods import highest_expressed_genes
from spapros.selection.selection_methods import random_selection
from spapros.selection.selection_methods import select_DE_genes
from spapros.selection.selection_methods import select_highly_variable_features
from spapros.selection.selection_methods import select_pca_genes
from spapros.selection.selection_methods import spca_feature_selection
from spapros.selection.selection_procedure import ProbesetSelector
from spapros.util.util import preprocess_adata


"""

Check below at __name__ == "__main__" for the description of this script


Notes:

We distinguish three types of parameters to specify selection run configurations:
- specific_params: They are method specific and keyword arguments for the according selection function
                   (defaults are given by function definitions of each method)
- general_params: They are method unspecific e.g. dataset, number of selected genes, ...
                  (defaults are given by `GENERAL_PARAMS` below)
- scheme_params: specific and general params always have one value refering to one configuration. However, scheme
                  params ease the way of defining multiple values for multiple configurations. They only occur in
                  experiment_configs but not in configs for each selection. E.g. `n_seeds` for random selections
                  creates a number of seeds for the method specific parameter `seed`. Potential `"seed"` values in
                  specific_params are overwritten.
                  Scheme params can also affect general params. E.g. `n_fraction_seeds`.
                  Even though default values are given as examples in SCHEME_PARAMS the params are ignored if not
                  set in the config file.
                  Scheme params are handled in `get_params_from_scheme_params`
"""

        
def spapros_selection_wrapper(adata,n,**kwargs):
    """Select a probeset with the spapros procedure
    A few default values are changed.
    """
    if "n_min_markers" not in kwargs:
        kwargs["n_min_markers"] = 0
    if "genes_key" not in kwargs:
        kwargs["genes_key"] = None
    selector = ProbesetSelector(adata,n=n,**kwargs,verbosity=0,n_jobs=-1)
    selector.select_probeset()
    return selector.probeset


METHODS = {
    "spapros": spapros_selection_wrapper,
    "pca": partial(select_pca_genes, inplace=False),
    "spca": partial(spca_feature_selection, inplace=False),
    "DE": partial(select_DE_genes, inplace=False),
    "random": partial(random_selection, inplace=False),
    "highest_expr": partial(highest_expressed_genes, inplace=False),
    "hvg": partial(select_highly_variable_features, inplace=False),
}

SCHEME_PARAMS = {
    "n_seeds": 1, # this is method specific (so far only method "random" supports this, some others should too...)
    "n_fraction_seeds": 1, # method unspecific, it affects the general param "fraction_seed" 
}

GENERAL_PARAMS = {
    "n": 100,
    "penalty_keys": [],
    "dataset": "small_data_raw_counts.h5ad",
    "data_path": "../package_dev/data/",
    "process_adata": ["norm","log1p"],
    "gene_subset": None,
    "obs_subset": None,
    "obs_fraction": None,  # NOT IMPLEMENTED, For measuring method stability
    "fraction_seed": None,  # NOT IMPLEMENTED, (obs_fraction and fraction_seed only work in combination)
    # for the fractions it might be more interesting to use disjoint subsets of the data!
    # on the other hand methods and statiscial power break with two few samples! (think you can better assess method
    # stability with the fraction+seed because you run that more often)
}

def get_params_from_scheme_params(scheme_params):
    """Create dictionary of general and specific params inferred from scheme params
    
    Arguments
    ---------
    scheme_params: dict
    
    Returns
    -------
    dict
    
    """
    params = {}
    for key, value in scheme_params.items():
        if value:
            if (key == "n_seeds"):
                np.random.seed(0)
                params["seed"] = list(np.random.choice(max([10000,value]), value, replace=False).tolist())
            elif (key == "n_fraction_seeds"):
                np.random.seed(0)
                params["fraction_seed"] = list(np.random.choice(max([10000,value]), value, replace=False).tolist())
            else:
                raise ValueError(f"Unknown scheme parameter: {key}")
    return params

def experiment_config_to_selection_configs(config_file,out_dir="./selection_configs"):
    """Convert config file for multiple selections to config files for each selection
    
    Config files for each selection are saved in `out_dir` as selection_{i}.yaml.
    Each selection_{i}.yaml has general and method specific parameters. The organisation of the `config_file` is more 
    complicated:
    You can define general and scheme parameters that apply to every method and then method specific parameters. 
    Additionally to method specific parameters you can again provide scheme and/or general parameters that are only 
    used for the given selection method and eventually overwrites configurations of the generally set parameters:
        scheme: (priority S2, G2)
            ...
        general: (priority G4)
            ...
        methods:
            pca:
                specific: (priority S3)
                    ...
                (scheme:) (priority S1, G1)
                (general:) (priority G3)
            ...
    
    Arguments
    ---------
    config_file: str
        Path to config yaml file for multiple selections.
    out_dir: str
        Directory where config for each selection is saved.
    """
    
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    
    # Load full config
    with open(config_file, "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    
    count = 0
    
    # Produce selection configs for each method separately
    for method in config["methods"]:
        
        spec_info = (config["methods"][method] is not None)
        
        ##############################################
        # Prepare general_params and specific_params #
        ##############################################
        
        # Get parameters based on scheme parameters        
        scheme_based_params = {}
        if ("scheme" in config) and config["scheme"]:
            scheme_based_params = get_params_from_scheme_params(config["scheme"])
        if spec_info and ("scheme" in config["methods"][method]) and config["methods"][method]["scheme"]:
            spec_scheme_based_params = get_params_from_scheme_params(config["methods"][method]["scheme"])
            # Add to or overwrite general scheme parameters with method specific scheme parameters
            scheme_based_params.update(spec_scheme_based_params)
        
        # Get general params. In case of doubled info the hierarchy is: 
        # default < general < method specific general < scheme based
        general_params = {key:[val] for key,val in GENERAL_PARAMS.items()}
        if ("general" in config) and config["general"]:
            for key, vals in config["general"].items():
                if (vals is not None) and (len(vals) > 0):
                    general_params[key] = vals
        if spec_info and ("general" in config["methods"][method]) and config["methods"][method]["general"]:
            for key, vals in config["methods"][method]["general"].items():
                if (vals is not None) and (len(vals) > 0):
                    general_params[key] = vals            
        for key,vals in scheme_based_params.items():
            if key in general_params:
                general_params[key] = vals
                
        # Get method specific params (scheme based params might overwrite defined params)
        specific_params = {}
        if spec_info and ("specific" in config["methods"][method]):
            for key,vals in config["methods"][method]["specific"].items():
                if (vals is not None) and (len(vals) > 0):
                    specific_params[key] = vals
        for key,vals in scheme_based_params.items():
            if key not in general_params:
                specific_params[key] = vals
        
        ##################################################################
        # Generate all possible configuration combinations and save them #
        ##################################################################
        cartesian_product = list(
            itertools.product(*[param_list for _, param_list in general_params.items()])  # type: ignore
        )
        general_configs = [{key: val for key, val in zip(general_params, val_list)} for val_list in cartesian_product]
        cartesian_product = list(
            itertools.product(*[param_list for _, param_list in specific_params.items()])  # type: ignore
        )
        method_configs = [{key: val for key, val in zip(specific_params, val_list)} for val_list in cartesian_product]        
        
        for g_config in general_configs:
            for m_config in method_configs:
                selection_config = {"method":method,"general":g_config,"specific":m_config}
                with open(os.path.join(out_dir,f"selection_{count}.yaml"), "w") as file:
                    yaml.dump(selection_config, file, allow_unicode=True)
                count += 1

                
def select(config_file,out_dir="./selections"):
    """Select gene sets according defined selection configuration
    
    Arguments
    ---------    
    config_file: str
        Path to config yaml file for a single selection.
    out_dir: str
        Directory where selection results and selection info are saved as csvs.
        
    """
    
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    
    # Load selection config
    with open(config_file, "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    general_params = config["general"]
    kwargs = config["specific"]
    method = config["method"]
    name = config_file.split("/")[-1].split(".")[0]
    
    print("Load data", flush=True)
    # Prepare data
    adata = sc.read(general_params["data_path"] + general_params["dataset"])
    print(adata, flush=True)
    preprocess_adata(adata, options=general_params["process_adata"], inplace=True)
    if general_params["gene_subset"] is not None:
        adata = adata[:,adata.var[general_params["gene_subset"]]]
    if ("obs_subset" in general_params) and (general_params["obs_subset"] is not None):
        adata = adata[adata.obs[general_params["obs_subset"]]]        
    print(adata, flush=True)
    
    # Initialize output files
    kwargs_keys = [k for k,v in kwargs.items() if not isinstance(v,dict)]
    kwdictargs_keys = [f"{k1}_{k2}" for k1,sub_dict in kwargs.items() if isinstance(sub_dict,dict) for k2 in sub_dict]
    param_keys = list(general_params) + kwargs_keys + kwdictargs_keys #list(kwargs) #
    param_keys.remove("process_adata")
    df_info = pd.DataFrame(
        columns=["method", "normalised", "log1p", "scaled", "time_seconds"] + param_keys
    )
    df_info.index.name = "set_id"
    df_sets = pd.DataFrame(index=adata.var.index)
    
    # Gene set selection
    print(f"Start selection for method {method}", flush=True)
    print(kwargs, flush=True)    
    start = timer()
    generated_selection = METHODS[method](adata, n=general_params["n"], **kwargs)  # type: ignore
    computation_time = timer() - start

    # Save method and computation time
    df_info.loc[name, ["method", "time_seconds"]] = [method,computation_time]
    
    # Save general params
    g_config_cols = [k for k in general_params if k in df_info.columns]
    g_config_values = [
        v if (not isinstance(v, list)) else "-".join(v) for k, v in general_params.items() if k in df_info.columns
    ]
    df_info.loc[name, g_config_cols] = g_config_values
    tmp = general_params["process_adata"]
    pp_options = [("norm" in tmp), ("log1p" in tmp), ("scale" in tmp)]
    df_info.loc[name, ["normalised", "log1p", "scaled"]] = pp_options
    
    # Save method specific params
    #kwarg_cols = [k for k in kwargs if k in df_info.columns]
    kwargs_values = [kwargs[k] if (not isinstance(kwargs[k], list)) else "-".join(kwargs[k]) for k in kwargs_keys]
    #kwarg_values = [v if (not isinstance(v, list)) else "-".join(v) for k, v in kwargs.items()]
    #df_info.loc[name, kwarg_cols] = kwarg_values
    df_info.loc[name, kwargs_keys] = kwargs_values
    kwdictargs_values = [kwargs[k1][k2] for k1,sub_d in kwargs.items() if isinstance(sub_d,dict) for k2 in sub_d]
    df_info.loc[name, kwdictargs_keys] = kwdictargs_values
    
    
    # Save selected genes
    df_sets[name] = generated_selection["selection"]

    # Save info and selection to files
    df_info.to_csv(os.path.join(out_dir,f"info_{name}.csv"))
    df_sets.to_csv(os.path.join(out_dir,f"{name}.csv"))
    # Additionally we save the selection specific gene infos (e.g. scores) in case we need to check them later
    generated_selection.to_csv(os.path.join(out_dir,f"scores_{name}.csv"))


def summarize_selections(selections_dir="./selections",out_dir="./"):
    """Combine multiple selection csvs to one
    
    Arguments
    ---------    
    selections_dir: str
        Directory where csvs of individual selections and infos of selections are saved.
    out_dir: str
        Directory where pooled are saved.
        
    """
    
    all_files = [f for f in os.listdir(selections_dir) if (f[-4:] == ".csv")]
    info_files = [os.path.join(selections_dir,f) for f in all_files if (f[:4] == "info")]
    files = [os.path.join(selections_dir,f) for f in all_files if (f[:4] != "info") & (f[:6] != "scores")]
    
    df = pd.concat([pd.read_csv(f,index_col=0) for f in files],axis=1)
    df_info = pd.concat([pd.read_csv(f,index_col=0) for f in info_files])
    
    # In case selections were on different datasets fill nan values
    df = df.fillna(False)
    
    df.to_csv(os.path.join(out_dir,"selections.csv"))
    df_info.to_csv(os.path.join(out_dir,"selections_info.csv"))


if __name__ == "__main__":
    """
    
    This script allows one to run a full selection in a step wise manner for parallelised workflows.
    
    step 0:
        Create multiple selection config files from one experiment config file.
    step 1:
         Run the selection for a given selection config file.
    step 2:
        Summarize the results of multiple runs of step 2 in one output (2 csv files)
    
    Example:
        # Create config files "./tmp_configs/selection_{i}.yaml"
        python select_script.py 0 ./selections_config.yaml ./tmp_configs
        
        # Run selections 
        python selection.py 1 ./tmp_configs/selections_0.yaml ./tmp_selections
        python selection.py 1 ./tmp_configs/selections_1.yaml ./tmp_selections
        python selection.py 1 ./tmp_configs/selections_2.yaml ./tmp_selections
        ...
        
        # Summarize results
        python selection.py 2 ./tmp_selections ./
    
    """
    
    step = int(sys.argv[1])
        
    print(f"Start selection step {step}", flush=True)
    if (step == 0):
        all_selections_config = sys.argv[2]
        configs_out_dir = sys.argv[3]
        experiment_config_to_selection_configs(
            config_file=all_selections_config,
            out_dir=configs_out_dir
        )
    elif (step == 1):
        selection_config = sys.argv[2]
        selections_out_dir = sys.argv[3]
        select(
            config_file=selection_config,
            out_dir=selections_out_dir
        )
    elif (step == 2):
        selections_dir = sys.argv[2]
        summary_out_dir = sys.argv[3]
        summarize_selections(
            selections_dir=selections_dir,
            out_dir=summary_out_dir
        )
    

