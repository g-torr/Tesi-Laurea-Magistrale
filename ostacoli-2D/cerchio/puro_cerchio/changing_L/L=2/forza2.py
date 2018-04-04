def forza2(start):
        fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))
	x,y=transpose(loadtxt("./posizione/dati_"+str(start)))
	N=len(fx)
	x=x[:10000]
	y=y[:10000]
	fx=fx[:10000]
	fy=fy[:100000]

        '''for i in arange(1,nf):
                temp_fx,temp_fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_x,temp_y=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx=append(fx,temp_fx)
		fy=append(fx,temp_fy)
		x=append(x,temp_x)
		y=append(x,temp_y)'''
	t=linspace(min(y),max(y),100)
	Fy=[]#chiamo Fy la forza cumulativa
	somma=0
	for k in arange (1,len(t)):
		for j in arange (1,len(x)):
			if ((fy[j]<t[k])&(fy[j]>t[k-1])):		
				somma+=fy[j]
		Fy=append(Fy,somma)
	return Fy	
