def stampa():
	lista=[]
	f=np.load("./R=0.5/F_per_part.npy")
	lista+=[f]
	f=np.load("./R=0.75/F_per_part.npy")
        lista+=[f]
	f=np.load("./R=1/F_per_part.npy")
	lista+=[f]
	f=np.load("./R=1.5/F_per_part.npy")
	lista+=[f]
	f=np.load("./R=2/F_per_part.npy")
	lista+=[f]
	f=np.load("./R=3/F_per_part.npy")
	lista+=[f]
	f=np.load("./R=4/F_per_part.npy")
	lista+=[f]
	f=np.load("./R=6/F_per_part.npy")
	lista+=[f]
	f=np.load("./R=8/F_per_part.npy")
	lista+=[f]
	'''R,p=np.load("./R=12/F_per_part.npy")
	lista+=[f]
	'''
	'''R,p=np.load("./R=16/pressione.npy")
	lista+=[f]
	R,p=np.load("./R=32/pressione.npy")
	lista+=[f]'''


	R=[0.5,0.75,1,1.5,2,3,4,6,8]
	return R,lista

