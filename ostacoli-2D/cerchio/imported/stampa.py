lista=[]
Raggi=[1,2,4,6,8,12,18,24,36,48,96,192]
for R in Raggi:
	R,p,var_p=np.load("./R="+str(R)+"/pressione.npy")
  	lista+=[R,p,var_p]
R,p,var_p=transpose(reshape(lista,[12,3]))
def stampa():
  figure(figsize=(4,3))
  cond=(1./R<0.1)
  plot(1./R[cond],p[cond]-1.,'o',label="Simulation",ms=6,mec="k",mew=1,c="w")
  fit=polyfit(1/R[cond],p[cond]-1., 1)
  plot(1./R[cond],polyval(fit,1/R[cond]),label="Fit",alpha=0.7)
  m=linspace(log10(1/max(R[cond])),log10(1/min(R[cond])));m=10**m
  plot(m,m,"r--",label="UCNA")
  fill_between(1/R[cond],p[cond]-1-var_p[cond],p[cond]-1+var_p[cond],label="error",color="gray",alpha=0.5)
  xlabel ("$\\mathscr{L}/R$",fontsize=10)
  ylabel("$\\left(\\hat{p}-\\hat{p}_0\\right)/\\hat{p}_0}$",fontsize=10)
  loglog()
  ylim(0.004,0.3)
  xlim(0.0045,0.1)
  tight_layout()
  xlabel ("$\\mathscr{L}/R$",fontsize=18)
  ylabel("$\\left(\\hat{p}-\\hat{p}_0\\right)/\\hat{p}_0}$",fontsize=18)
  return R,p,var_p

def stampa_all():
  figure(figsize=(4,3))
  cond=(1./R<0.1)
  cond2=(1/R<2)
  #R=R[cond2];p=p[cond2];var_p=var_p[cond2];
  plot(1./R[cond2],p[cond2]-1.,'o',label="Simulation",ms=6,mec="k",mew=1,c="w")
  fit=polyfit(1/R[cond],p[cond]-1., 1)
  plot(1./R[cond],polyval(fit,1/R[cond]),label="Fit",alpha=0.7)
  m=linspace(log10(1/max(R)),log10(1/min(R)));m=10**m
  plot(m,m,"r--",label="UCNA")
  fill_between(1/R[cond2],p[cond2]-1-var_p[cond2],p[cond2]-1+var_p[cond2],label="error",color="gray",alpha=0.5)
  xlabel ("$\\mathscr{L}/R$",fontsize=10)
  ylabel("$\\left(\\hat{p}-\\hat{p}_0\\right)/\\hat{p}_0}$",fontsize=10)
  loglog()
  ylim(0.004,1)
  xlim(0.0045,2)
  tight_layout()
  xlabel ("$\\mathscr{L}/R$",fontsize=18)
  ylabel("$\\left(\\hat{p}-\\hat{p}_0\\right)/\\hat{p}_0}$",fontsize=18)


  return R,p,var_p

