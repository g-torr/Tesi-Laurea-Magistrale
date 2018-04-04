import scipy.special
r1=3.
r2=1.
def stampa():
	D,f_1,p_bordo_1=transpose(loadtxt("./tau=0.03/sim_vs_D.dat"))
	D,f_2,p_bordo_2=transpose(loadtxt("./tau=0.06/sim_vs_D.dat"))
	D,f_3,p_bordo_3=transpose(loadtxt("./tau=0.12/sim_vs_D.dat"))
	
	t=linspace(-4,3,100)
	d=10**t
	#subplot(311)
	figure(figsize=[4,3])
	'''plot(t,pi/(2*pi*sqrt(t)*(1+sqrt(t))+pi),color="b",label="UCNA,$D=1$")
	plot(t,pi/(2*pi*sqrt(4*t)*(1+sqrt(4*t))+pi),color="magenta",label="UCNA,$D=4$")'''
	plot(d,pi*3/(normalization(3,1,sqrt(d*0.03))),color="b")
	plot(d,pi*3/(normalization(3,1,sqrt(d*0.06))),color="magenta")
	plot(d,pi*3/(normalization(3,1,sqrt(d*0.12))),color="gray")
	plot(D,(1-p_bordo_1),"s",label="$\\tau=0.03$",ms=5,mec="b",mew=2,c="w",alpha=0.7)	
	plot(D,(1-p_bordo_2),"o",label="$\\tau=0.06$",ms=5,mec="m",mew=2,c="w",alpha=0.7)
	plot(D,(1-p_bordo_3),"^",label="$\\tau=0.12$",ms=5,mec="gray",mew=2,c="w",alpha=0.7)
#	axvline(0.5,color="red",ms=3)
	legend(loc="best",numpoints=1,frameon=False)
	xlabel("$D$",fontsize=8)	
	ylabel("$\\rho_{Bulk}$",fontsize=8)
	semilogx()
	tight_layout()
	xlabel("$D$",fontsize=18)	
	ylabel("$\\rho_{Bulk}$",fontsize=18)


	l=linspace(-2,1,100)
	l=10**l
	figure(figsize=[4,3])
	#subplot(312)
	plot(l,pi*r1*r2/(normalization(r1,r2,l)),color="k",label="SCA")
	plot(sqrt(D*0.03),(1-p_bordo_1),"s",label="$\\tau=0.03$",ms=5,mec="b",mew=2,c="w",alpha=0.7)
	plot(sqrt(D*0.06),(1-p_bordo_2),"o",label="$\\tau=0.06$",ms=5,mec="m",mew=2,c="w",alpha=0.7)
	plot(sqrt(D*0.12),(1-p_bordo_3),"^",label="$\\tau=0.12$",ms=5,mec="gray",mew=2,c="w",alpha=0.7)
	semilogx()
	legend(loc="best",numpoints=1,frameon=False)
	xlabel("$\\mathscr{L}$",fontsize=8)	
	ylabel("$\\rho_{Bulk}$",fontsize=8)
	tight_layout()
	xlabel("$\\mathscr{L}$",fontsize=18)	
	ylabel("$\\rho_{Bulk}$",fontsize=18)

	#subplot(313)



	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	figure(figsize=[4,3])
	
	
	plot(sqrt(D*0.03),f_1/arc_lenght,"s",label=" $\\tau=0.03$",ms=5,mec="b",mew=2,c="w",alpha=0.7)
	plot(sqrt(D*0.06),f_2/arc_lenght,"o",label=" $\\tau=0.06$",ms=5,mec="m",mew=2,c="w",alpha=0.7)
	plot(sqrt(D*0.12),f_3/arc_lenght,"^",label=" $\\tau=0.12$",ms=5,mec="gray",mew=2,c="w",alpha=0.7)
	plot(l,l**2/0.03*(arc_lenght+2*pi*l)/(normalization(r1,r2,l)*arc_lenght),color="b")
	plot(l,l**2/0.06*(arc_lenght+2*pi*l)/(normalization(r1,r2,l)*arc_lenght),color="m")
	plot(l,l**2/0.12*(arc_lenght+2*pi*l)/(normalization(r1,r2,l)*arc_lenght),color="gray")
	xlim(0.015,10)
	ylim(0.0005,10)

	xlabel("$\\mathscr{L}$",fontsize=8)	
	ylabel("$p$",fontsize=8)
	loglog()
	legend(loc="best",numpoints=1,frameon=False)
	tight_layout()
	xlabel("$\\mathscr{L}$",fontsize=18)	
	ylabel("$p$",fontsize=18)
	xlim(0.01,11)
	ylim(0.001,20)	
	legend(loc=(0.,0.55),numpoints=1,frameon=False)
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	figure(figsize=[4,3])
	
	
	plot(sqrt(D*0.03),f_1/D/(arc_lenght/(pi*r1*r2)),"s",label=" $\\tau=0.03$",ms=5,mec="b",mew=2,c="w",alpha=0.7)
	plot(sqrt(D*0.06),f_2/D/(arc_lenght/(pi*r1*r2)),"o",label=" $\\tau=0.06$",ms=5,mec="m",mew=2,c="w",alpha=0.7)
	plot(sqrt(D*0.12),f_3/D/(arc_lenght/(pi*r1*r2)),"^",label=" $\\tau=0.12$",ms=5,mec="gray",mew=2,c="w",alpha=0.7)
	plot(l,(arc_lenght+2*pi*l)/(normalization(r1,r2,l)*arc_lenght/(pi*r1*r2)),color="b")
	plot(l,(arc_lenght+2*pi*l)/(normalization(r1,r2,l)*arc_lenght/(pi*r1*r2)),color="m")
	plot(l,(arc_lenght+2*pi*l)/(normalization(r1,r2,l)*arc_lenght/(pi*r1*r2)),color="gray")
	
	xlabel("$\\mathscr{L}$",fontsize=8)	
	ylabel("$p/p_0$",fontsize=8)
	legend(loc="best",numpoints=1,frameon=False)
	semilogx()
	tight_layout()
	xlabel("$\\mathscr{L}$",fontsize=18)	
	ylabel("$p/p_0$",fontsize=18)

def normalization(r1,r2,L):
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	return 	pi*r1*r2+L*arc_lenght+L**2*2*pi	

				
