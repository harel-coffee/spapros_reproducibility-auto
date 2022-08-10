
1. Create Conda Environments:
        
The initial environment: 
      
    conda env create --prefix ./environments/venv_yaml --file environments/venv_yaml.yaml

For each method: 

      conda env create --prefix ./environments/venv_method --file environments/venv_method.yaml

Note: for scmer and nsforest use venv_gene_selection.yaml      


2. Charliecloud container for R methods
    
    TODO: Adjust script selection_script_charliecloud.sh to be independent of ICB cluster


3. How to run external selections:

- Set up 1. and 2. 
- Adjust all configs in config_yamls_<...> to the correct data path
- run bash run_ext_selections.sh