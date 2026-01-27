# Uncertainty Quantification for Estimates of Extreme Isolines of Bivariate Distributions

Project to explore estimation and uncertainty quantification of extreme isolines of bivariate distributions, with application to bivariate climate events.

## How to..

### Clone this Repository

1. Run `git clone https://github.com/jbbutler/isolines_uq.git` in desired directory

### Install Conda Environment and Expose to JupyterHub

1. If at NERSC, run `module load conda`
1. In the project home directory, run `conda env create -f environment.yml`
1. Activate the environment with `conda activate isolinesR`
1. Start an R session with `R`
1. Run the following commands line-by-line:
```
> ename <- Sys.getenv('CONDA_DEFAULT_ENV')
> dname <- trimws(paste("R", getRversion(), Sys.getenv("CONDA_PROMPT_MODIFIER")))
> IRkernel::installspec(name=ename, displayname=dname)
> quit()
```
You should now see this environment as available to use in a Jupyter notebook.

### Run Tutorial Notebook

1. Navigate to `/isolines_uq/notebooks/tutorials/tube_construction.ipynb`. Open the notebook.
1. Make sure the notebook is running the `R 4.2.3 (isolinesR)` kernel (upper righthand corner of notebook). If not, go to `Kernel > Change Kernel` and find it in the dropdown.