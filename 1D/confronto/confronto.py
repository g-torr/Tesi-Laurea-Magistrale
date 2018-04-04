import os
from  pylab import *
D=0.1
tau=1.
size=8.6
L=sqrt(D*tau)
def confronto():	
	figure(figsize=(4,3))
	x_comp,m_comp,err=np.load("./compenetranti/forza_tot.npy")
	er_comp=abs(m_comp-m_comp[-1])/m_comp[-1]
	x_cnew=x_comp[:7];x_cnew=append(x_cnew,x_comp[7:19:2]);x_cnew=append(x_cnew,x_comp[19:])
	er_cnew=er_comp[:7];er_cnew=append(er_cnew,er_comp[7:19:2]);er_cnew=append(er_cnew,er_comp[19:])
	plot(x_cnew,er_cnew,"o",ms=6,mec="black",mew=2,c="w",label="Hard Potential",alpha=0.6)
	fill_between(x_comp,er_comp - 8e-06/m_comp[-1],er_comp+ 8e-06/m_comp[-1],color=(.4,.4,.4),alpha=0.5)


	x_soff,m_soff=np.load("./soffice/forza_tot.npy")
	er_soff=abs(m_soff-m_soff[-1])/m_soff[-1]
	x_cnew=x_soff[:7];x_cnew=append(x_cnew,x_soff[7:19:2]);x_cnew=append(x_cnew,x_soff[19:])
	er_cnew=er_soff[:7];er_cnew=append(er_cnew,er_soff[7:19:2]);er_cnew=append(er_cnew,er_soff[19:])
	plot(x_cnew,er_cnew,"o",ms=6,mec="b",mew=2,c="w",label="Soft Potential",alpha=0.6)
	fill_between(x_soff,er_soff - 8e-06/m_soff[-1],er_soff+ 8e-06/m_soff[-1],color=(.4,.4,.9),alpha=0.5)
	xlabel("$\\Delta t /\\tau$",fontsize=14)
	ylabel("Relative variation",fontsize=14)
	semilogx()
	xlim(0,0.3)
	ylim(-0.001)

def confronto2():
	figure(figsize=(4,3))
	x_comp,m_comp,err=np.load("./compenetranti/forza_tot.npy")
	x_soff,m_soff=np.load("./soffice/forza_tot.npy")
	a,b=np.load("/home/giuseppe/Documents/myCUDA/1D/compenetranti/compenetranti_dt/forza0.800000.npy")
	x_comp=insert(x_comp,0,a);m_comp=insert(m_comp,0,b);
	plot(x_comp,m_comp,lw=2,c="m",label="Hard Potential")
	plot(x_soff,m_soff,c="c",lw=2,label="Soft Potential")
	t=linspace(-4.,0.,100)
	t=10**t
	plot(t,2*D/(2*L+size)*1**t,"--",label="UCNA")
	fill_between([0.25,0.8], 0.020,0.026,color="gray",alpha=0.5)
	xlabel("$\\Delta t/\\tau$",fontsize=12)
	ylabel("$\\langle |\Phi '| \\rangle$",fontsize=12)
	semilogx()
	xlim(0,0.8)
	ylim(0.020,0.026)

	
