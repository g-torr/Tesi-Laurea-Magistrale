	
def forza2(start,stop):
	area,arco=np.load("area.npy")
        fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))
	x,y=transpose(loadtxt("./posizione/dati_"+str(start)))
	N=len(fy)
	fy=fy/N
	for i in arange(start+1,stop+1):
                temp_fx,temp_fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_x,temp_y=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx=append(fx,temp_fx/N)
		x=append(x,temp_x)
		fy=append(fy,temp_fy/N)
		y=append(y,temp_y)
	ratio=float(len(fy[(fy!=0)|(fx!=0)]))/float(len(fy)) #particelle sul bordo/tutte le particelle
	h,b=histogram(y,bins=arange(min(y)-0.02,max(y)+0.02,0.00001),weights=fy)
	h=area*h/((1-ratio)*(stop-start+1))#normalizzo rispetto alla densita' di bulk
	sh=cumsum(h)
	b=(b[:-1]+b[1:])/2
	figure(0)
	title ("$F_y(y)$")
	xlabel("$y$")
	ylabel ("$F_y$")
	np.save("F.npy",array([b,h]))
	plot(b,h)
	figure(1)
	title ("cumulative $F_y(y)$")
	xlabel("$y$")
	ylabel ("$cum F_y$")
	plot(b,sh)
	np.save("F_cumulative.npy",array([b,sh]))

