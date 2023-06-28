# Graded Sensitivity to Structure and Meaning throughout the Human Language Network

This repository provides scripts to support reproduction of the results reported in
Shain, Kean, et al. (2023).

## Installation

Python dependencies can be managed through the Anaconda framework.
Download and install conda (https://www.anaconda.com/), then initialize a new environment called `nlength`:

    conda env create -f conda.yml

Before running any scripts, activate the `nlength` environment:

    conda activate nlength

To run the statistical analyses, you will also need to install `R` (https://www.r-project.org/) and add the `lme4` library:

    R
    >>> install.packages('lme4')

## Data

Data are publicly available for download at https://osf.io/fduve/.



## Usage

This repository provides scripts that must be run sequentially to produce the results reported in the paper.
All scripts must be run from the root directory of this repository.

### 1. Initialize the path to the data downloaded from OSF
    python -m nlength.set_data_path

### 2. Compute the contrast estimates

    python -m nlength.contrasts

This will create a subdirectory called `contrasts` with tables containing all contrasts needed for analyses.
It will also create a subdirectory called `plots` with visualizations of these contrasts.

### 3. (Optional) Plot the key contrasts by region of interest

    python -m nlength.plot

This will add bar plots to the `plots` directory.

### 4. Fit group-level models

    nlength/regress_l2.R

This will create a subdirectory called `glm` with LM/E models and summary reports for each test in this study.

### 5. Run group-level tests

    nlength/test_l2.R

This will create a subdirectory called `lrt` with summary reports of all likelihood ratio tests.

### 6. Tabulate test results

    python -m nlength.signif_table

This will create a table called `signif.csv` with FDR-correct p-values.

### 7. (Optional) Compute the supplementary SWJN analyses

    python -m nlength.swjn

This will create a table called `swjn.csv` and add `SWJN.png` to the `plots` directory.

### 8. (Optional) Report the results of the language localizer validation analyses

    python -m nlength.validate_langloc

### 9. (Optional) Generate brain surface plots

Download the Surf Ice software (https://www.nitrc.org/projects/surfice/) to this directory and run `nlength/surfice.py` from the Surf Ic GUI.

## References

Shain, Kean, et al. (2023). Graded Sensitivity to Structure and Meaning throughout the Human Language Network. _bioRxiv_.