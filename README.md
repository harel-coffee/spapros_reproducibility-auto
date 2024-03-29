# spapros_reproducibility

This repository contains the code to reproduce the analyses and figures of the manuscript "Probeset selection for targeted spatial transcriptomics". The manuscript is currently under review. Preprint at [bioRxiv](https://doi.org/10.1101/2022.08.16.504115).

The repository is structured as follows:
- preprint folder: contains the scripts and notebooks to reproduce the analyses and figures till the preprint version. This includes Figures 1c, 2c,d,e 3, S1, S2a, S3, S4, S5, S6, S7, S9b,c.
- revision folder: contains the scripts and notebooks to reproduce the analyses and figures for the revised version of the manuscript. This includes Figures 4, S2b, S8, S9a, S10, S11, S12, S13.
- For the revised version we ran extensive benchmarks with our snakemake pipeline (release tag 0.1.0). The configurations for rerunning the experiments can be found [there](https://github.com/theislab/spapros-smk) at `spapros-smk/configs/` with the file names `config_rev.yaml`, `config_rev_bm.yaml`, `config_rev_bm_probes.yaml`.

Please do not hesitate to open an issue in case of any questions or issues.


