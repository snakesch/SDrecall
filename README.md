# SDrecall

## TODO
- [ ] Three separate directories at top level (preparation, phasing, shell_utils)
- [x] Need a way to locate target shell scripts (Alternatively: source shell_utils.sh; check_bam_validity?)
- [ ] **Preparation done. Moving `sd_pairs.py` and `convert_nodes_into_hierachical_beds.py` soon.**
- [ ] Numba not yet in env
- [ ] README not fixed
- [ ] conda YAML lacking some libraries (e.g. intervaltree)
- [ ] Extra documentation for ad-hoc tools for intrinsic variant calling
- [ ] Review license
- [ ] Cleanup code


## Installation
### Using conda/mamba
Users should first clone this repository to a local directory.

For mamba users, create an environment from YAML.
```{bash}
mamba env create -f ./env/SDrecall.yml
mamba activate SDrecall
```

## Usage
(help message of driver code)

### Inputs
Users need to provide a QC-ed BAM file and an optional BED regions of selected variant recall regions.

### Outputs

### Using docker/singularity
Given the long list of dependencies of SDrecall, we are still working on a docker file / singularity recipe. Any contributions are most welcome. 

## Code development and feature request
SDrecall is under active development. We welcome all kinds of suggestions and collaborations.

## Contact and correspondance
Xingtian Yang (yangyxt@hku.hk), Louis She (snakesch@connect.hku.hk)
