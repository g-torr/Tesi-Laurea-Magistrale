import os
import sys
import re

for dirname, dirnames, filenames in os.walk('.'):
    break
'''dati=[]
for name in dirnames:
	os.chdir(name)
	l=float(re.sub("L=","",name))
	execfile("analisi.py")
	dati+=[l,bulk_probability(9,20)]
	os.chdir("../")
	#dati+=[l,p]
L=array(dati[::2])
rh=dati[1::2]


np.save("bulk_density.npy",[L,rh])
'''
def density():
	L,rh=np.load("bulk_density.npy")
	l=linspace(-2,2,100)
	l=10**l
	plot(l,1/(pi*(1+l*(2*l+sqrt(2*pi)))))
	plot(L,rh,"o")
def density_distribution(start,stop):
	x,y=transpose(loadtxt("./posizione/dati_"+str(start)))
	for i in arange(start+1,stop+1):
		temp_x,temp_y=transpose(loadtxt("./posizione/dati_"+str(i)))
		x=append(x,temp_x)
		y=append(y,temp_y)
	r=sqrt(x**2+y**2)
	t=linspace(-3.5,0,400)
	t=(1.-10**(t))[::-1]
	h,b=histogram(r, bins=t, normed=False)
	bins=(b[:-1]+b[1:])/2.
	plot(1-bins, h/(2*pi*bins*diff(b)*len(r)))#note that normed=False is important, in this way h is equal to the count the number of particle in the bin, then I divide by the surface, i.e. 2*pi*bins, and the number of particles, i.e. len(r). Binning is not perfect, so I may miss some particle in the binning, but doesn't matter in this way and the result is consistent for the density of bulk.

def plot_4_distributions():
	Ls=[2,0.5,0.1,0.025]
	f, axarr = plt.subplots(2, 2)
	start=10
	stop=19
	rows,column=indices((2,2))
	rows=rows.flatten()
	column=column.flatten()
	for l in Ls:
		x,y=transpose(loadtxt("L="+str(l)+"/posizione/dati_"+str(start)))
		for i in arange(start+1,stop+1):
			temp_x,temp_y=transpose(loadtxt("L="+str(l)+"/posizione/dati_"+str(i)))
			x=append(x,temp_x)
			y=append(y,temp_y)
		r=sqrt(x**2+y**2)
		t=linspace(-3.5,0,400)
		t=(1.-10**(t))[::-1]
		h,b=histogram(r, bins=t, normed=False)
		bins=(b[:-1]+b[1:])/2.
		i=Ls.index(l)
		axarr[rows[i],column[i]].plot(1-bins, h/(2*pi*bins*diff(b)*len(r)),label="L="+str(l))
		axarr[rows[i],column[i]].legend(loc="best")
		if l<1.:
			axarr[rows[i],column[i]].axvline(l,ls="--",color="r")
		axarr[1,0].set_xlabel("r")
		axarr[1,1].set_xlabel("r")
		axarr[0,0].set_ylabel("Density")
		axarr[1,0].set_ylabel("Density")
def distribution():
	L=[]
	item=array([(float(name[2:]),name) for name in dirnames],dtype=[("L",np.float),("name","S7")])
	for l,name in sort(item,order="L"):
		os.chdir(name)
		l=float(re.sub("L=","",name))
		density_distribution(10,15)
		L+=[l]
		os.chdir("../")
	legend(L,loc="best")
def pressure():
	dati=[]
	for name in dirnames:
		L=float(re.sub("L=","",name))
		try:
			R,p,var_p=np.load("./"+str(name)+"/pressione.npy")

		except IOError:
			continue
		dati+=[L,p]
	l=linspace(0,1,100)
	figure(figsize=(4,3))
	plot(l, sqrt(pi/2)*l,"--",label="Prediction")
	L=array(dati[::2])
	press=array(dati[1::2])
	plot(L,press-1,'o',label="Simulation",ms=6,mec="m",mew=2,c="w")
	loglog()
	xlabel("$\mathscr{L}/R$",fontsize=10)
	ylabel("$\\left(\\hat{p}-\\hat{p}_0\\right)/\\hat{p}_0}$",fontsize=10)
	tight_layout()
	xlabel("$\mathscr{L}/R$",fontsize=15)
	ylabel("$\\left(\\hat{p}-\\hat{p}_0\\right)/\\hat{p}_0}$",fontsize=15)
	legend(loc="best",numpoints=1,frameon=False)
	ylim(0,1.7)
	return L,press #return L and pressure
