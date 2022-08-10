import sys
import pandas as pd
import scanpy as sc
import yaml
from pathlib import Path
from spapros.evaluation import ProbesetEvaluator
from spapros.util.util import preprocess_adata

"""

Python script for step-wise evaluation of probesets. The separation of steps enables parallelised evaluation of multiple sets.

Arguments (see below at __name__ == "__main__")
---------
results_dir:
    Location to save evaluation results.
evaluation_config:
    Path to configuration yaml file.
probeset_csv:
    Path to csv file that contains the probesets.
metrics:
    Metrics to evaluate.
steps:
    The evaluation steps (can be a subset of [0,1,2,3]).
set_ids:
    Set ids from ``probeset_csv`` to evaluate.


"""



def get_genes(set_id, var_names, probesets_file="./probesets.csv"):
    """ """
    #selection = pd.read_csv(probesets_file, usecols=["index", set_id], index_col=0)
    selection = pd.read_csv(probesets_file, index_col=0)
    genes = [g for g in selection.loc[selection[set_id]].index.to_list() if g in var_names]
    return genes

def evaluate(
    results_dir="./probeset_evaluation",    
    evaluation_config="./evaluation_config.yaml",
    probeset_csv="./probesets.csv",    
    metrics=["cluster_similarity"],
    steps=[0,1,2,3],    
    set_ids=None,
):
    """Run evaluation for defined configs, probesets and evaluation steps 
    
    results_dir: str
        Save evaluation results here. The reference results are saved in `results_dir+"reference".
    metrics_config: str
        Path to yaml file with metric configurations.
    data_config: str
        Path to yaml file with data set configurations.
    probeset_csv: str
        Path to csv file with probe sets.
    metrics: list of strs
        Metrics to evaluate.
    steps: list of ints
        Steps that are computed for the given metrics. Note when not starting with 0 then result files of previous 
        stepts might be needed in `results_dir`. 
    set_ids: list of strs
        List of probesets in `probeset_csv` that are evaluated.
        
    """
    
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    
    # Load configs
    with open(evaluation_config, "r") as file:
        params = yaml.load(file, Loader=yaml.FullLoader)
    data_params = params["data"]
    metric_params = params["metrics"]
    
    # Load and process data
    adata = sc.read(data_params["data_path"] + data_params["dataset"])
    preprocess_adata(adata, options=data_params["process_adata"])
    if ("gene_subset" in data_params) and (data_params["gene_subset"] is not None):
        adata = adata[:,adata.var[data_params["gene_subset"]]]
    if ("obs_subset" in data_params) and (data_params["obs_subset"] is not None):
        adata = adata[adata.obs[data_params["obs_subset"]]]
    print(adata, flush=True)    
    
    # Load probe sets
    gene_sets = {set_id:get_genes(set_id, adata.var_names, probesets_file=probeset_csv) for set_id in set_ids}
    
    # Evaluation    
    evaluator = ProbesetEvaluator(adata,
                                  scheme="custom",
                                  results_dir=results_dir,
                                  celltype_key=data_params["celltype_key"] if ("celltype_key" in data_params) else None,
                                  marker_list=metric_params["marker_corr"]["marker_list"] if ("marker_corr" in data_params) else None,
                                  metrics_params=metric_params,
                                  metrics=metrics,
                                  reference_name=data_params["name"], 
                                  reference_dir=None,  
                                  verbosity=0,
                                  n_jobs=-1,
                                 )    
    
    if 0 in steps:
        evaluator.compute_or_load_shared_results()
    if 1 in steps:
        for set_id,genes in gene_sets.items():
            evaluator.evaluate_probeset(genes, set_id=set_id, pre_only=True)
    if 2 in steps:
        for set_id,genes in gene_sets.items():
            evaluator.evaluate_probeset(genes, set_id=set_id, update_summary=False)
    if 3 in steps:
        evaluator.summary_statistics([set_id for set_id in gene_sets])


if __name__ == "__main__":
    """
    
    Example of running this script:
    
    python eval_script.py ./probeset_evaluation ./evaluation_config.yaml ./probesets.csv 
    marker_corr-gene_corr 0123 set1 set2 set3 ... setn
    
    """
     
    results_dir = sys.argv[1]
    evaluation_config = sys.argv[2]
    probeset_csv = sys.argv[3]
    metrics = sys.argv[4].split("-")
    steps = [int(s) for s in list(sys.argv[5])]
    set_ids = sys.argv[6:]
    
    print(f"{results_dir=}", flush=True)
    print(f"{evaluation_config=}", flush=True)
    print(f"{probeset_csv=}", flush=True)
    print(f"{metrics=}", flush=True)
    print(f"{steps=}", flush=True)
    print(f"{set_ids=}", flush=True)
    
    
    evaluate(
        results_dir=results_dir,
        evaluation_config=evaluation_config,
        probeset_csv=probeset_csv,
        metrics=metrics,
        steps=steps,
        set_ids=set_ids,
    )
    