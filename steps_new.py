from optparse import OptionParser
import math
import numpy as np
import pandas as pd
import scipy.stats
from pysnptools.snpreader import Bed

# step 2
def generate_betas(M,h1,h2,rho,p00,p10,p01,p11): 
	
	# first generate causal effect sizes gamma as outlined in UNITY paper
	
	# covariance matrix of gammas
        cov1 = h1 / (M * (p11 + p10))  
        cov2 = math.sqrt(h1)*math.sqrt(h2)*rho / (M * p11)  
        cov3 = cov2
        cov4 = h2 / (M * (p11 + p01))
        cov = [[cov1, cov2], [cov3, cov4]]

	# generating gammas from bivariate normal distribution 
        print('Simulating gammas...')
        gammas = np.random.multivariate_normal([0,0], cov, M)

	# now generate causal indicator vector for each trait

	# determine which traits each variant is causal for by drawing from multinomial distribution  
        print('Simulating Cs...')
        Cs = np.random.multinomial(1, [p00,p10,p01,p11], M)
	
	# causal indicator vectors
        C1 = np.zeros(M)
        C2 = np.zeros(M)

	# inferring causal indicator vectors from draws of multinomial distribution 
        print('Getting causal vectors...')
        C1[np.where(Cs[:,1])] = 1
        C1[np.where(Cs[:,3])] = 1
        C2[np.where(Cs[:,2])] = 1
        C2[np.where(Cs[:,3])] = 1
	
	# true effect size beta is product of gamma and C 
        print('Getting betas...')
        betas1 = gammas[:,0] * C1
        betas2 = gammas[:,1] * C2

        return betas1, betas2

# step 3
def get_beta_tildes(bed_file,mean_std_file,betas1,betas2,h1,h2,chunk_size_snp):

	# reading bed file 
	G = Bed(bed_file,count_A1=False)
	# reading file with means, standard deviation for each SNP 
	mean_std = pd.read_csv(mean_std_file,delimiter='\t')

	# dimensions of genotype matrix
	N = G.row_count # number of individuals 
	M = G.col_count # number of SNPs

	# dot products of standardized matrix and betas
	GB1 = np.zeros(N)
	GB2 = np.zeros(N)

	# standardizing genotype matrix and taking dot product with betas (chunk_size_snp at a time)
	for i in range(0,M,chunk_size_snp):
		
		# standardizing 
		G_sub = G[:, i : (i + chunk_size_snp)].read().val # current chunk	
		mean_sub = mean_std['mean'][i : i + chunk_size_snp].values # means of SNPs corresponding to current chunk
		std_sub = mean_std['std'][i: (i + chunk_size_snp)].values # standard deviations of SNPs corresponding to current chunk 
		nanidx = np.where(np.isnan(G_sub)) # finding NaNs in genotype matrix 
		G_sub[nanidx] = mean_sub[nanidx[1]] # setting NaNs to mean of corresponding SNP
		G_sub_std = np.nan_to_num((G_sub - mean_sub) / std_sub) # standardizing chunk 
		
		# dot product
		betas1_sub = betas1[i : (i + chunk_size_snp)] # trait 1 effect sizes of SNPs corresponding to current chunk
		betas2_sub = betas2[i : (i + chunk_size_snp)] # trait 2 effect sizes of SNPs corresponding to current chunk 
		GB1 += np.dot(G_sub_std,betas1_sub) # dot product for trait 1
		GB2 += np.dot(G_sub_std,betas2_sub) # dot product for trait 2 

	# re-scaling to have variance of dot product equal to heritability 
	var_GB1 = np.var(GB1)  
	var_GB2 = np.var(GB2) 
	k1 = h1 / var_GB1
	k2 = h2 / var_GB2
	beta_tildes1 = math.sqrt(k1) * betas1 # re-scaled effect sizes for trait 1 
	beta_tildes2 = math.sqrt(k2) * betas2 # re-scaled effect sizes for trait 2	

	return beta_tildes1, beta_tildes2 

# steps 4-5
def get_epsilon_tildes(N,h1,h2):
	
	# generating environmental noise as outlined in UNITY paper 
	e1 = np.random.normal(0,math.sqrt(1-h1),N) # noise for trait 1 
	e2 = np.random.normal(0,math.sqrt(1-h2),N) # noise for trait 2 
	
	# re-scaling to have variance of noise equal to 1 - heritability 
	var_e1 = np.var(e1)
	var_e2 = np.var(e2)
	l1 = (1 - h1) / var_e1
	l2 = (1 - h2) / var_e2
	e_tildes1 = math.sqrt(l1) * e1 # re-scaled noise for trait 1
	e_tildes2 = math.sqrt(l2) * e2 # re-scaled noise for trait 2

	return e_tildes1, e_tildes2 
	
# step 6
def get_beta_hats(bed_file,mean_std_file,beta_tildes1,beta_tildes2,e_tildes1,e_tildes2,chunk_size_snp,chunk_size_ind):
	
	# reading bed file
	G = Bed(bed_file,count_A1=False)
	# reading file with means, standard deviation for each SNP 
	mean_std = pd.read_csv(mean_std_file,delimiter='\t')

	# dimensions of genotype matrix 
	N = G.row_count # number of individuals 
	M = G.col_count # number of SNPs
	
	# phenotype vectors (y = G * beta tilde + epsilon)
	# initializing with noise; dot product will be added in loop
	y1 = e_tildes1 
	y2 = e_tildes2

	# standardizing genotype matrix (alternatively, could output it in step 3 and input it here) and taking dot product with beta tildes  
	for i in range(0,M,chunk_size_snp):
	
		# standardizing 
		G_sub = G[:, i : (i + chunk_size_snp)].read().val # current chunk
		mean_sub = mean_std['mean'][i : i + chunk_size_snp].values # means of SNPs corresponding to current chunk
		std_sub = mean_std['std'][i: (i + chunk_size_snp)].values # standard deviations of SNPs corresponding to current chunk
		nanidx = np.where(np.isnan(G_sub)) # finding NaNs in genotype matrix  
		G_sub[nanidx] = mean_sub[nanidx[1]] # setting NaNs to mean of corresponding SNP 
		G_sub_std = np.nan_to_num((G_sub - mean_sub) / std_sub) # standardized chunk 

		# dot product
		beta_tildes1_sub = beta_tildes1[i : (i + chunk_size_snp)] # trait 1 effect sizes of SNPs corresponding to current chunk
		beta_tildes2_sub = beta_tildes2[i : (i + chunk_size_snp)] # trait 2 effect sizes of SNPs corresponding to current chunk
		y1 += np.dot(G_sub_std,beta_tildes1_sub) # dot product for trait 1  
		y2 += np.dot(G_sub_std,beta_tildes2_sub) # dot product for trait 2
		
	# "OLS" estimates for beta (beta hat)
	beta_hats1 = np.zeros(M)
	beta_hats2 = np.zeros(M)

	# standardizing genotype matrix (again) and finding "OLS" solution for betas (beta hat = (tranpose(G) * y) / N)
	for i in range(0,N,chunk_size_ind):
		
		# chunks are now "horizontal": all SNPs are used, only select number of individuals used 	
	
		# standardizing 
		G_sub = G[i : (i + chunk_size_ind), :].read().val # current chunk 
		means = mean_std['mean'].values # means of all SNPs
		std = mean_std['std'].values # standard deviations of all SNPs
		nanidx = np.where(np.isnan(G_sub)) # finding NaNs in genotype matrix 
		G_sub[nanidx] = means[nanidx[1]] # setting NaNs to mean of correspoinding SNP
		G_sub_std = np.nan_to_num((G_sub - means) / std) # standardizing chunk 

		# transpose(G) * y
		y1_sub = y1[i : (i + chunk_size_ind)] # trait 1 phenotypes of individuals corresponding to current chunk 
		y2_sub = y2[i : (i + chunk_size_ind)] # trait 2 phenotypes of individuals corresponding to current chunk 
		beta_hats1 += np.dot(G_sub_std.transpose(),y1_sub) # dot product for trait 1 
		beta_hats2 += np.dot(G_sub_std.transpose(),y2_sub) # dot product for trait 2 

	# dividing by N 
	beta_hats1 = beta_hats1 / N 
	beta_hats2 = beta_hats2 / N 

	return beta_hats1, beta_hats2 

def main():

	parser = OptionParser()
	parser.add_option("--h1",dest="h1") # heritability of trait 1 
	parser.add_option("--h2",dest="h2") # heritability of trait 2 
	parser.add_option("--rho",dest="rho") # genetic correlation 
	parser.add_option("--p00",dest="p00") # proportion of variants that are not causal for either trait
	parser.add_option("--p10",dest="p10") # proportion of variants that are causal for only trait 1
	parser.add_option("--p01",dest="p01") # proportion of variants that are causal for only trait 2
	parser.add_option("--p11",dest="p11") # proportion of variants that are causal for both traits
	parser.add_option("--seed",dest="seed") 
	parser.add_option("--sim_num",dest="sim_num",default=None) # only necessary for naming output file 
	parser.add_option("--chr_i",dest="chr_i",default=None) # only necessary if data is of only one chromosome; otherwise, assumes all data is provided 
	parser.add_option("--bed_file",dest="bed_file") # file must end in ".bed"
	parser.add_option("--mean_std_file",dest="mean_std_file") # tab-delimited file with columns named "mean", "std" 
	parser.add_option("--chunk_size_snp",dest="chunk_size_snp",default=1000) # number of SNPs to use in one computation (when using all individuals)  
	parser.add_option("--chunk_size_ind",dest="chunk_size_ind",default=600) # number of individuals to use in one computation (when using all SNPs)

	(options, args) = parser.parse_args()	

	h1 = float(options.h1)
	h2 = float(options.h2)
	rho = float(options.rho)
	p00 = float(options.p00)
	p10 = float(options.p10)
	p01 = float(options.p01)
	p11 = float(options.p11)
	seed = int(options.seed)
	if options.sim_num is None:
        	sim_num = seed # if no simulation number is provided, seed number is used to name output file 
	else:
		sim_num = int(options.sim_num)
	chr_i = int(options.chr_i)
	bed_file = options.bed_file
	mean_std_file = options.mean_std_file
	chunk_size_snp = int(options.chunk_size_snp)
	chunk_size_ind = int(options.chunk_size_ind)

	np.random.seed(seed)
		
	# dimensions of genotype matrix
	N = Bed(bed_file,count_A1=False).row_count # number of individuals
	M = Bed(bed_file,count_A1=False).col_count # number of SNPs	

	# random effect sizes (beta)
	print('Generaitng betas...')
	(betas1, betas2) = generate_betas(M,h1,h2,rho,p00,p10,p01,p11)

	# random effect sizes standardized to have variance equal to heritability (beta tilde)
	print('Getting beta tildes...')
	(beta_tildes1, beta_tildes2) = get_beta_tildes(bed_file,mean_std_file,betas1,betas2,h1,h2,chunk_size_snp)

	# random environmental noise standardized to have variance equal to 1 - heritability (epsilon tilde)
	print('Generating epsilons and getting epsilon tildes...')
	(e_tildes1, e_tildes2) = get_epsilon_tildes(N,h1,h2)

	# OLS estimates for beta 
	print('Getting beta hats...')
	(beta_hats1, beta_hats2) = get_beta_hats(bed_file,mean_std_file,beta_tildes1,beta_tildes2,e_tildes1,e_tildes2,chunk_size_snp,chunk_size_ind)

	# reading bim file (alternatively, could turn this into an input; make it up to the user to provide correct bim file)
	df = pd.read_csv('/u/scratch/r/ruthjohn/mario_ruthie/ukbb_sumstats/meta_data/filter4.bim',delimiter='\t',names=['CHR','SNP','CM','BP','A1','A2'])
	
	if chr_i is None:	
		# if no chromosome is specified, assumes data from all chromosomes is provided  
		df1 = df.copy()
		df2 = df.copy()
	else:
		# if chromosome is specified, filters bim file by given chromosome 
		df1 = df[df.CHR==chr_i]
		df2 = df[df.CHR==chr_i]

	# column with sample sizes (N)
	N_col = np.ones(M)*N

	# column with z-scores (Z)
	z1 = beta_hats1 / math.sqrt(1/N)
	z2 = beta_hats2 / math.sqrt(1/N)
	
	# column with p-values (P)
	p1 = (1 - scipy.stats.norm.cdf(abs(z1)))*2
	p2 = (1 - scipy.stats.norm.cdf(abs(z2)))*2

	# adding columns to dataframes 
	df1['N'] = N_col
	df1['Z'] = z1
	df1['BETA_STD'] = beta_hats1
	df1['P'] = p1

	df2['N'] = N_col
	df2['Z'] = z2
	df2['BETA_STD'] = beta_hats2
	df2['P'] = p2

	# writing files 
	out_file1 = '/u/scratch/r/ruthjohn/mario_ruthie/unity_v2.0/results_new_sims/sims_files/p_{}_h_{}_{}_rho_{}_{}.gwas1'.format(p11,h1,h2,rho,sim_num)
	out_file2 = '/u/scratch/r/ruthjohn/mario_ruthie/unity_v2.0/results_new_sims/sims_files/p_{}_h_{}_{}_rho_{}_{}.gwas2'.format(p11,h1,h2,rho,sim_num)	
	df1.to_csv(out_file1,sep= ' ',index=False)
	df2.to_csv(out_file2,sep=' ',index=False)
	
if __name__ == "__main__":
	main()
