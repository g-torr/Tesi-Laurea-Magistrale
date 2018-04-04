from pylab import*
print("it is specifically designed for quenced disorder") 
def forza2(start,stop):
	area,arco=transpose(loadtxt("area.txt"))
        fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))
	x,y=transpose(loadtxt("./posizione/dati_"+str(start)))
	N=len(fy)
	fy=fy
	for i in arange(start+1,stop+1):
                temp_fx,temp_fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_x,temp_y=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx=append(fx,temp_fx)
		x=append(x,temp_x)
		fy=append(fy,temp_fy)
		y=append(y,temp_y)
	r0=0.7
	F=sqrt(fx**2+fy**2)
	h,b=histogram(y,bins=arange(-1.,max(y)+0.02,0.0005),weights=fy) #it is the yforce density distribution, density is normalized at N
	h=h/N/(stop-start+1)#normalizzo rispetto al numero di particelle
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

def termalizzazione(index):
	area,arco=transpose(loadtxt("area.txt"))
        fx,fy=transpose(loadtxt("./forza/forza_"+str(index)))
	x,y=transpose(loadtxt("./posizione/dati_"+str(index)))
	N=len(fy)
	r0=0.7
	F=sqrt(fx**2+fy**2)
	h,b=histogram(y,bins=arange(-1.,max(y)+0.02,0.0005),weights=fy) #it is the yforce density distribution, density is normalized at N
	h=h/N #normalizzo rispetto al numero di particelle
	sh=cumsum(h)
	b=(b[:-1]+b[1:])/2
	figure(0)
	title ("$F_y(y)$")
	xlabel("$y$")
	ylabel ("$F_y$")
	np.save("F"+str(index)+".npy",array([b,h]))
	plot(b,h)
	figure(1)
	title ("cumulative $F_y(y)$")
	xlabel("$y$")
	ylabel ("$cum F_y$")
	plot(b,sh)
	np.save("F_cumulative"+str(index)+".npy",array([b,sh]))

