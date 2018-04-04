def genera(R=3,r=0.5,theta=pi/8,n=50):
	t=linspace(theta,pi-theta,n)      #arco piccolo
	T=linspace((pi-theta),2*pi+theta,n*int(R/r))	#arco grande
	box=(R-r)/sin(theta)/2.
	x=r*cos(t)
	y=r*sin(t)+(R-r)/sin(theta)-box
	X=R*cos(T)
	Y=R*sin(T)-box
	plot(x,y)
	plot(X,Y)
	mx=1.1*max(abs(concatenate((x,y,X,Y))))
	#xlim(-mx,mx)
	#ylim(-mx,mx)
	ax = gca()
	ax.set_aspect(1)
	x1=linspace(-5,5,100)
	y1=R*sin(theta)-1./tan(theta)*x1+ R*(cos(theta))**2/sin(theta) -box
	y2=R*sin(theta)+1./tan(theta)*(x1+ R*cos(theta)) -box	
	plot(x1,y1)
	plot(x1,y2)
	xtot=concatenate((x,X))
	ytot=concatenate((y,Y))
	savetxt("input.dat",transpose(array([xtot,ytot])).astype(float32))
	plot(xtot,ytot,lw=2)
	print("il numero di punti  generati is "+str(len(xtot)))
