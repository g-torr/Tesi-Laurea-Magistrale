def stampa():
	lista=[]
	"""R,p=np.load("./R=0.5/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=1/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=1.5/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=2/pressione.npy")
	lista+=[R,p]"""
	R,p=np.load("./R=3/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=4/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=6/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=8/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=12/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=16/pressione.npy")
	lista+=[R,p]

	r=linspace(0,0.5,100)
	plot(r,sqrt(pi/2)*r,"--",label="Prediction")
	R,p=transpose(reshape(lista,[6,2]))
	plot(1/R,p-1,'o',ms=10,mec="k",mew=3,c="w")
	fit=polyfit(1/R,p-1, 1)
	plot(1/R,polyval(fit,1/R),label="fit")
	xlabel ("$\\frac{1}{R}$")
	ylabel("Pressure")
	return R,p

