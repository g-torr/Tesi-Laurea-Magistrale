def monta():
	dt=[0.00326]
	for i in arange(0,5):
	    dt+=[dt[-1]/2.]
	lista=[]
	dt=array(dt)
	for i in arange(0,len(dt)): 
		dts= '{0:.6f}'.format(dt[i])
		nome="forza"+dts+".npy"
		#nome_err="errore"+dts+".npy"
		x,m=np.load(nome)
		#x,err=np.load(nome_err)
		#lista+=[[x,m,err]]
		lista+=[[x,m]]
	#lista=array(lista)
	#return(lista)	
	x,m=transpose(lista)
	#fill_between(x,m-err,m+err,color=(.7,.7,.7))
	plot(x,m,"-o")
	return x,m
