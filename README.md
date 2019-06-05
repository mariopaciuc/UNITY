# UNITY Genotype Data Simulations

# Getting Started

The required versions of various packages can be found in ```pip_freeze.txt```. 

Run ```https://github.com/mariopaciuc/Unity_simulations.git``` to donlowad the Unity simulator and ```cd Unity_simulations``` to proceed. 

# Simulating Data

Data is simulated with ```steps_new.py```. To simulate data, the user must specify the proportion heritability of each trait (h1, h2), the geentic correlation between the traits (rho), and the proportion of causal SNPs for each trait (p00, p10, p01, p11). The user must also provide plink-style data in ```.bed/.bim/.fam``` format and a tab-delimited file with the mean and standard deviation of each SNP with column names "mean" and "std". 

The user also has the option of specifying chunk sizes to split up the compuations in a more memory-efficient manner. One chunk size relates to the number of SNPs and another to the number of individuals. The remaining optional inputs are seed, simulation number, and chromosome number. If running multiple simulations with the same set of parameters, the user should specify the seed (and change it for each simulation)--otherwise, all the simulations will yield the same data. The simulation number is required only for naming the output file--if no number is provided, the seed is used to differentiate among different simulations. 

The simulated data can be found in ```results_new_sims/sims_flies``` in a directory according to its true value of p11, rho, and h1 and h2. The simulated effect sizes are under the column "BETA_STD" as part of a dataframe with other useful information. 

Example:
```python steps_new.py --h1 .25 --h2 .25 --rho 0 --p00 .97 --p10 .01 --p01 .01 --p11 .01 --seed 1 --sim_num 1 --chunk_size_snp 1000 --chunk_size_ind 600 --bed_file chr3.bed --bim_file chr3.bim --mean_std_file mean_std_file.txt```

