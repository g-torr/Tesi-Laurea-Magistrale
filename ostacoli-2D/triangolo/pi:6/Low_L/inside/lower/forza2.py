from pylab import*
L=sqrt(0.1*0.06)
theta=pi/6
def forza2(start,stop):
        fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))
	x,y=transpose(loadtxt("./posizione/dati_"+str(start)))
	N=len(fy)
	for i in arange(start+1,stop+1):
                temp_fx,temp_fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_x,temp_y=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx=append(fx,temp_fx)
		x=append(x,temp_x)
		fy=append(fy,temp_fy)
		y=append(y,temp_y)
	
	'''cond=y<(cos(theta)/sin(theta)*(x-L)+1/sin(theta)-1)
	cond2=y< (-cos(theta)/sin(theta)*(x+L)+1/sin(theta)-1)
	cond3=y> (-sin(theta)+1)/cos(theta)*(x+L)
	cond4=y> (sin(theta)-1)/cos(theta)*(x-L)'''
	
	F=sqrt(fx**2+fy**2)
	tx=2*rand(1e5)-1
	ty=2*rand(1e5)-1
	area=count_nonzero(filtering_function(tx,ty))/float(len(tx))*4.
	rho=count_nonzero(filtering_function(x,y))/area #densita' di bulk
	h,b=histogram(y,bins=arange(-1.,max(y)+0.005,0.0005),weights=fy) #it is the yforce density distribution, density is normalized at N
	h=h/rho#normalizzo rispetto alla densita' di bulk
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
	r=sqrt(x**2+y**2)
	F=sqrt(fx**2+fy**2)
	rho=count_nonzero(r<r0)/(pi*r0**2) #densita' di bulk
	h,b=histogram(y,bins=arange(-1.,max(y)+0.02,0.0005),weights=fy) #it is the yforce density distribution, density is normalized at N
	h=h/rho#normalizzo rispetto alla densita' di bulk
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

def filtering_function(x,y):
	cond=y<(cos(theta)/sin(theta)*(x-2*L)+1/sin(theta)-1)
	cond2=y< (-cos(theta)/sin(theta)*(x+2*L)+1/sin(theta)-1)
	cond3=y> (-sin(theta)+1)/cos(theta)*(x+2*L)
	cond4=y> (sin(theta)-1)/cos(theta)*(x-2*L)
	return (cond & cond2 &(cond3 | cond4))

