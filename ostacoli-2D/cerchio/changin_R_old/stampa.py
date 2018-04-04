def stampa():
	lista=[]
	R,p=np.load("./R=0.5/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=0.75/pressione.npy")
        lista+=[R,p]
	R,p=np.load("./R=1/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=1.5/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=2/pressione.npy")
	lista+=[R,p]
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
	'''R,p=np.load("./R=16/pressione.npy")
	lista+=[R,p]
	R,p=np.load("./R=32/pressione.npy")
	lista+=[R,p]'''


	R,p=transpose(reshape(lista,[10,2]))
	return R,p

