def genera(n):#inserire il numero di lati del poligono
	t=linspace(0,2*pi,n+1)#[:-1]
	x=cos(t)
	y=sin(t)
	plot(x,y)
	x=x[:-1]
	y=y[:-1]	
	savetxt("input2.dat",transpose(array([x,y])).astype(float32))
	
	
