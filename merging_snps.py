# removes SNPs from file based on LD pruning 
# input: dataset obtained from simulations (sim_file), list of SNPs to keep (snps_file)
# output: dataset with only desired SNPs in same file as one input 

from optparse import OptionParser
import pandas as pd

parser = OptionParser()
parser.add_option("--sim_file",dest="sim_file")
parser.add_option("--snps_file",dest="snps_file")

(options, args) = parser.parse_args()

sim_file = options.sim_file
snps_file = options.snps_file

sims_df = pd.read_csv(sim_file,delimiter=' ')
snps_df = pd.read_csv(snps_file,names=['SNP'])
pruned_df = sims_df.merge(snps_df,how='inner',on='SNP')

# overwrites original sims file
out_file = sim_file 
pruned_df.to_csv(out_file,sep=' ',index=False)
