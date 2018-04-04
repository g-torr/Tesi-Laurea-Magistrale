def genera(r1=1,r2=1,n=100):#inserire il numero di lati del poligono
	t=linspace(0,2*pi,n+1)#[:-1]
	x=r1*cos(t)
	y=r2*sin(t)
	plot(x,y)
	x=x[:-1]
	y=y[:-1]	
	savetxt("input.dat",transpose(array([x,y])).astype(float32))
	k=r1*r2/(((r1**2-r2**2)*sin(t)**2+r2**2)**1.5)
	np.save("curvatura",k)
