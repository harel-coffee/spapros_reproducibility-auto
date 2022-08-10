import os
import yaml
import sys
import pandas as pd
import scanpy as sc


###########################
# Main selection function #
###########################

def run_selection(config_file, out_dir):
    """Run a selection according the params in the config_file yaml

    Arguments
    ---------
    config_file: str
        Path to config yaml.
    out_dir: str
        Path to directory where selected gene set and selection info are saved.

    Save
    ----
    1. os.path.join(out_dir,f"info_{name}.csv")
        Table with infos of selection parameters
    2. os.path.join(out_dir,f"{name}.csv")
        Table with
            index: adata.var_names
            column "selection" (bool): Selected genes
    """

    # Load selection config
    with open(config_file, "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    # Define parameters
    general_params = config["general"]
    kwargs = config["specific"]
    if kwargs is None:
        kwargs = {}
    method = config["method"]
    # name = config_file.split("/")[-1].split(".")[0]  # not working on windows
    name = os.path.basename(config_file).split(".")[0]
    n = general_params["n"]

    # Load data
    adata = sc.read(os.path.join(general_params["data_path"],general_params["dataset"]))
    if ("gene_subset" in general_params) and (general_params["gene_subset"] is not None):
        adata = adata[:,adata.var[general_params["gene_subset"]]].copy()
    if ("obs_subset" in general_params) and (general_params["obs_subset"] is not None):
        adata = adata[adata.obs[general_params["obs_subset"]]].copy()

    # Run selection
    # SCGENEFIT
    if method == "scgenefit":
        from selection_methods.gene_selection_scgenefit import select_genes_scgenefit
        selection, computation_time = select_genes_scgenefit(n,adata=adata,**kwargs)

    # NSFOREST
    elif method == "nsforest":
        from selection_methods.gene_selection_nsforest import gene_selection_nsforest
        # labl = kwargs["label"]
        # n_clust = adata.obs[labl].cat.categories.shape[0]
        selection, computation_time = gene_selection_nsforest(n, adata, **kwargs)  # int(np.ceil(n / n_clust))

    # SCPNMF
    elif method == "scpnmf":
        from selection_methods.gene_selection_scpnmf import select_genes_scpnmf
        tmp_dir = os.path.join(out_dir, "tmp")
        conda_env = config["venv"]
        adata = os.path.join(general_params["data_path"],general_params["dataset"])  # need to delete the anndata because otherwise the h5ad file is locked and can't be read by the R script
        if not os.path.exists(tmp_dir):
            os.umask(0)
            os.makedirs(tmp_dir, 0o777)
            os.chmod(tmp_dir, 0o777)
        # r_exe = "/usr/bin/Rscript"
        selection, computation_time = select_genes_scpnmf(n, adata, output_path=os.path.join(tmp_dir, "scpnmf_tmp.tsv"), **kwargs, conda_env=conda_env)
        adata = sc.read(adata)

    # SCMER
    elif method == "scmer":
        from selection_methods.gene_selection_scmer import select_genes_scmer
        selection, computation_time = select_genes_scmer(n, adata, **kwargs)

    # SMASH
    elif method == "smash":
        from selection_methods.gene_selection_smash import select_genes_smash
        selection, computation_time = select_genes_smash(n, adata, **kwargs)
        if 'method' in kwargs:
            name = name + "_" + kwargs["method"]
        else:
            name = name + "_DNN"

    # ASFs
    elif method == "asfs":
        from selection_methods.gene_selection_asfs import select_genes_asfs
        selection, computation_time = select_genes_asfs(n, adata, **kwargs)

    # geneBasis
    elif method == "genebasis":
        from selection_methods.gene_selection_genebasis import select_genes_genebasis
        tmp_dir = os.path.join(out_dir, "tmp")
        conda_env = config["venv"]
        if not os.path.exists(tmp_dir):
            os.umask(0)
            os.makedirs(tmp_dir, 0o777)
            os.chmod(tmp_dir, 0o777)
        selection, computation_time = select_genes_genebasis(n=n,
                                                             adata=adata,
                                                             tmp_dir=tmp_dir,
                                                             conda_env=conda_env,
                                                             **kwargs)

    # selfE
    elif method == "selfe":
        from selection_methods.gene_selection_selfe import select_genes_selfe
        tmp_dir = os.path.join(out_dir, "tmp")
        conda_env = config["venv"]
        if not os.path.exists(tmp_dir):
            os.umask(0)
            os.makedirs(tmp_dir, 0o777)
            os.chmod(tmp_dir, 0o777)
        selection, computation_time = select_genes_selfe(n=n,
                                                         adata=adata,
                                                         tmp_dir=tmp_dir,
                                                         conda_env=conda_env,
                                                         **kwargs)

    # COSG
    elif method == "cosg":
        from selection_methods.gene_selection_cosg import select_genes_cosg
        selection, computation_time = select_genes_cosg(n, adata, **kwargs)

    # Triku
    elif method == "triku":
        from selection_methods.gene_selection_triku import select_genes_triku
        selection, computation_time = select_genes_triku(n, adata, **kwargs)



    ###############
    # Output csvs #
    ###############
    # The following code just saves parmeters to one csv and the selected genes to another csv
    # It looks a little bit complicated since arugments that are dictionaries are also supported (needed that for our selection procedure....)
    # I think for most methods that's not needed, so when creating the analogous R script that could be simpler. And even here you can simplify it if you want.

    # Initialize output files
    kwargs_keys = [k for k,v in kwargs.items() if not isinstance(v,dict)]
    kwdictargs_keys = [f"{k1}_{k2}" for k1,sub_dict in kwargs.items() if isinstance(sub_dict,dict) for k2 in sub_dict]
    param_keys = list(general_params) + kwargs_keys + kwdictargs_keys
    df_info = pd.DataFrame(
        columns=pd.Series(["method", "time_seconds"] + param_keys).unique()
    )
    df_info.index.name = "set_id"
    df_sets = pd.DataFrame(index=adata.var.index)

    # Save method and computation time
    df_info.loc[name, ["method", "time_seconds"]] = [method,computation_time]

    # Save general params
    g_config_cols = [k for k in general_params if k in df_info.columns]
    g_config_values = [
        v if (not isinstance(v, list)) else "-".join(v) for k, v in general_params.items() if k in df_info.columns
    ]
    df_info.loc[name, g_config_cols] = g_config_values

    # Save method specific params
    kwargs_values = [kwargs[k] if (not isinstance(kwargs[k], list)) else "-".join(kwargs[k]) for k in kwargs_keys]
    df_info.loc[name, kwargs_keys] = kwargs_values
    kwdictargs_values = [kwargs[k1][k2] for k1,sub_d in kwargs.items() if isinstance(sub_d,dict) for k2 in sub_d]
    df_info.loc[name, kwdictargs_keys] = kwdictargs_values


    # Save selected genes
    df_sets[name] = False
    df_sets[name].iloc[selection.idx] = True

    # Save info and selection to files
    df_info.to_csv(os.path.join(out_dir,f"info_{name}.csv"))
    df_sets.to_csv(os.path.join(out_dir,f"{name}.csv"))



if __name__ == "__main__":

    config_file = sys.argv[1]
    out_dir = sys.argv[2]

    run_selection(config_file,out_dir)

