import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")
import scipy.signal
import sys
import scipy
import scipy.io
def fx(x,m,c):
	return m*x+c
def gt(t,tau,beta):
	w=(t/tau)**beta
	return np.exp(-w)

for i2 in range(len(sys.argv)):
	if sys.argv[i2]=="-cn":
		chain_number=int(sys.argv[i2+1])
	elif sys.argv[i2]=="-an":
		atom_number=int(sys.argv[i2+1])
	elif sys.argv[i2]=="-t":
		t_step=float(sys.argv[i2+1])
	elif sys.argv[i2]=="-s":
		tot=int(sys.argv[i2+1])
	else:
		pass

atom_number = int(atom_number)
chain_number = int(chain_number)
ac = int(atom_number/chain_number)
sq_X=np.zeros((ac,chain_number))
for j in range(chain_number):
	sq_X_a=np.fromfile("mean_sq_V_"+str(j), dtype=float)
	for i in range(ac):
		sq_X[i][j]=sq_X_a[i]

sq_X_mean=np.zeros(ac)
for i in range(ac):
	for j in range(chain_number):
		sq_X_mean[i]=sq_X_mean[i]+sq_X[i][j]/chain_number
#tot=1519
#print(sq_X_mean)
file_name2="sq_V_mean.txt"
writefile = open(file_name2,'w')
for j in range(ac):
	writefile.write(str(sq_X_mean[j]))
	#print(sq_X_mean[j])
	writefile.write("\n")
corr_a=np.zeros((tot,ac))
for j in range(chain_number):	
	matrix=np.fromfile("rousev_1_"+str(j),dtype=float)
	count=0
	for i2 in range(ac):
		for k in range(tot):
			corr_a[k][i2]=corr_a[k][i2]+matrix[k+i2*tot]/chain_number

ln_corr_a=np.zeros((tot,ac))
file_name2="corr_V.txt"
writefile = open(file_name2,'w')
for i in range(tot):
	writefile.write(str(i*t_step))
	writefile.write("\t")
	for j in range(ac):
		writefile.write(str(corr_a[i][j]))
		writefile.write("\t")
		ln_corr_a[i][j]=np.log(corr_a[i][j])
	writefile.write("\n")

ms=['o','s','^','*','p','d','x','v','>','<','h']

plt.figure()
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.gcf().set_size_inches(4,3,forward=True)
plt.subplots_adjust(left=0.18,bottom=0.15)
for j in range(100):
	plt.plot(j*t_step,corr_a[j][1],'o',markerfacecolor='w',markeredgecolor='b',markersize=4)
plt.savefig("V1.png",dpi=300)







