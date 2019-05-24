from optparse import OptionParser
import numpy as np
import dask
import dask.array as da 
from pandas_plink import read_plink

# have to figure out parser with python3, for now will just write in inputs in script

parser = OptionParser()
(options, args) = parser.parse_args()
parser.add_option("--plink_file", dest="plink_file")
parser.add_option("--chunk_size", dest="chunk_size")
parser.add_option("--out_file", dest="out_file")

#plink_file = options.plink_file
plink_file='../../../sample_plink_files/chrX'
(bim, fam, G) = read_plink(plink_file)

N = G.shape[0] 
M = G.shape[1] 
#chunk_size = options.chunk_size
chunk_size = M
G = G.rechunk((N,chunk_size))

means_vect = np.zeros(ncol)
sd_vect = np.zeros(ncol)

for i in range(0,M,chunk_size):
	print('looping for i=')
	print(i)
	G_sub = G[:, i:(i + chunk_size)]
	print('getting means...')
	sub_mean = G_sub.mean(0)
	print('getting standard deviations...')
	sub_sd = G_sub.std(0)
	#nanidx = np.where(np.isnan(G_sub))
	#G_sub[nanidx] = sub_mean[nanidx[1]]
	print('adding means to vector...')
	means_vect[i:(i + chunk_size)] = sub_mean
	print('adding standard deviations to vector...')
	sd_vect[i:(i + chunk_size)] = sub_sd
	
print('Getting matrix...')	
X_std = (G.compute() - means_vect) / sd_vect # standardized matrix 

#out_file = options.out_file	
#np.save(out_file,X_std) 	   
