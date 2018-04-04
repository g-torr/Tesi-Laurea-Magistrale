def monta(dt):
	lista=[]
	dt=array(dt)
	for i in arange(0,len(dt)): 
		dts= '{0:.6f}'.format(dt[i])
		nome="forza"+dts+".npy"
		x,y=np.load(nome)
		lista+=[[x,y]]
	#lista=array(lista)
	#return(lista)	
	x,y=transpose(lista)
	plot(x,y,"-o")
	return x,y
