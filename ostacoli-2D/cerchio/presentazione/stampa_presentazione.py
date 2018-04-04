import itertools
marker	= itertools.cycle(('o', 's', '^')) 
color	= itertools.cycle(((0.8,0,0), (0.8,0,0.8), (0,0.8,0.8), (0.1,0.5,0.)))
labels	= itertools.cycle(("$\\mathscr{L}=1$","$R=1,D=1$","$R=1,D=4$"))
def stampa():
	lista=[];dati=[]
	#Raggi=[1,2,4,6,8,12,18,24,36,48,96,192]
	Raggi=[2,4,6,8,12,18,24,36,48,96,192]
	L=1.
	for R in Raggi:
		R,p,var_p=np.load("/home/giuseppe/Documents/myCUDA/ostacoli-2D/cerchio/imported/R="+str(R)+"/pressione.npy")
  		lista+=[L/R,p,var_p]
	L_R,p,var_p=transpose(reshape(lista,[11,3]))
	figure(figsize=(4,3))
	grafica(L_R,p/1.-1.,var_p)
	dati=lista
	lista=[]
	L=[0.1,0.2,0.4,0.5,0.7,0.05,0.25,0.025]
	for l in L:
		R,p,var_p=np.load("/home/giuseppe/Documents/myCUDA/ostacoli-2D/cerchio/poligoni/D=1/L="+str(l)+"/pressione.npy")
   		lista+=[l/R,p,var_p]
	L_R,p,var_p=transpose(reshape(lista,[8,3]))
	grafica(L_R,p/1.-1.,var_p)
	dati=append(dati,lista)
	lista=[]
	tau=[0.01,0.04,0.08,0.25,0.0025,0.16]
	for t in tau:
		R,p,var_p=np.load("/home/giuseppe/Documents/myCUDA/ostacoli-2D/cerchio/poligoni/D=4/tau="+str(t)+"/pressione.npy")
   		lista+=[sqrt(dot(t,4))/R,p/4.,var_p/4.]
	L_R,p,var_p=transpose(reshape(lista,[6,3]))
	grafica(L_R,p-1.,var_p)
	dati=append(dati,lista)
	m=linspace(0,1,1000)
	plot(m,m,"r--",label="SCA")
	xlim(0.004,1.2)
	ylim(0.005,1.2)
	L_R,p,var=transpose(reshape(dati,(25,3)))
	L_R,p=sort([L_R,p])
	cond=(L_R<0.2)
	'''fit=polyfit(L_R[cond],p[cond]-1., 1)
	plot(L_R,polyval(fit,L_R),"-.k",label="Fit",alpha=0.7)'''
	x = L_R[cond]
	x=x[:,np.newaxis]
	a, _, _, _ = np.linalg.lstsq(x,p[cond]-1.) 
	plot(x,a*x,"-.k",label="Fit",alpha=0.7)
	legend(loc=(0,0.38),numpoints=1,frameon=False,fontsize=12)
	axvline(0.2,ls="--",color="gray")

	return L_R,p
def grafica(L_R,p,var_p):
  mec=color.next()
  plot(L_R,p,marker.next(),label=labels.next(),ms=6,mec=mec,mew=1,c="w")
  '''fit=polyfit(L/R,p, 1)
  plot(L/R,polyval(fit,L/R),label="Fit",alpha=0.7)'''
  errorbar(L_R,p,var_p,ecolor=mec,fmt=None)
  xlabel ("$\\mathscr{L}/R$",fontsize=10)
  ylabel("$\\hat{p}/\\hat{p}_0}$",fontsize=10)
  tight_layout()
  xlabel ("$\\mathscr{L}/R$",fontsize=18)
  ylabel("$(\\hat{p}-\\hat{p}_0)/\\hat{p}_0}$",fontsize=18)
  return L_R,p,var_p

