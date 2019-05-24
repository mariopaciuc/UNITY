# simulating betas

from optparse import OptionParser
import math
import numpy as np
import dask
import dask.array as da

def generate_betas(M,h1,h2,rho,p00,p10,p01,p11):

	sigma_sq1 = 1-h1
	sigma_sq2 = 1-h2

	cov1 = sigma_sq1 / (M * (p11 + p10))
	cov2 = math.sqrt(sigma_sq1)*math.sqrt(sigma_sq2)*rho / (M * p11)
	cov3 = cov2
	cov4 = sigma_sq2 / (M * (p11 + p01))
	cov = [[cov1, cov2], [cov3, cov4]]

	print('Simulating gammas...')
	gammas = np.random.multivariate_normal([0,0], cov, M)
	gammas = da.from_array(gammas)
	print('Simulating Cs...')
	Cs = np.random.multinomial(1, [p00,p10,p01,p11], M)
	Cs = da.from_array(Cs)

	C1 = np.zeros(M)
	C2 = np.zeros(M)

	# vectors indicating causal SNPs
	print('Getting causal vectors...')
	C1[np.where(Cs[:,1])] = 1
	C1[np.where(Cs[:,3])] = 1
	C2[np.where(Cs[:,2])] = 1
	C2[np.where(Cs[:,3])] = 1

	print('Getting betas...')
	betas1 = gammas[:,0] * C1
	betas2 = gammas[:,1] * C2

	return betas1, betas2

def get_beta_tildes(X_std,betas1,betas2,h1,h2):
	
	XB1 = X_std.dot(betas1)
	XB2 = X_std.dot(betas2)
	var_XB1 = XB1.var()
	var_XB2 = XB2.var()
	k1 = h1 / var_XB1
	k2 = h2 / var_XB2

	beta_tildes1 = math.sqrt(k1) * betas1
	beta_tildes2 = math.sqrt(k2) * betas2
	beta_tildes1 = beta_tildes1.compute()
	beta_tildes2 = beta_tildes2.compute()
		
	return beta_tildes1, beta_tildes2

def get_epsilon_tildes(N,h1,h2):
	
	e1 = da.random.normal(0,math.sqrt(1-h1),N)
	e2 = da.random.normal(0,math.sqrt(1-h2),N)
	var_e1 = e1.var()
	var_e2 = e2.var()
	l1 = h1 / var_e1
	l2 = h2 / var_e2

	e_tildes1 = math.sqrt(l1) * e1
	e_tildes2 = math.sqrt(l2) * e2
	e_tildes1 = e_tildes1.compute()
	e_tildes2 = e_tildes2.compute()
	
	return e_tildes1, e_tildes2
	
def get_beta_hats(X_std,beta_tildes1,beta_tildes2,e_tildes1,e_tildes2):
	
	y1 = X_std.dot(beta_tildes1) + e_tildes1
	y2 = X_std.dot(beta_tildes2) + e_tildes2
		

parser = OptionParser()
(options, args) = parser.parse_args()
parser.add_option("--h1",dest="h1")
parser.add_option("--h2",dest="h2")
parser.add_option("--rho",dest="rho")
parser.add_option("--p00",dest="p00")
parser.add_option("--p10",dest="p10")
parser.add_option("--p01",dest="p01")
parser.add_option("--p11",dest="p11")
parser.add_option("--X_std_file",dest="X_std_file")

N = 300000
M = 500000
#h1 = options.h1
#h2 = options.h2
#rho = options.rho
#p00 = options.p00
#p10 = options.p10
#p01 = options.p01
#p11 = options.p11

h1=.5
h2=.5
rho=.01
p00=.93
p10=.01
p01=.01
p11=.05
