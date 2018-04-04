def genera(R=3,theta=pi/6,n=50):
#	t=linspace(theta,pi-theta,n)      #arco piccolo
	T=linspace((pi-theta),2*pi+theta,n)	#arco grande
#	box=(R-r)/sin(theta)/2.
	box=0
	x=0
	y=0-box
	X=R*cos(T)
	Y=R*sin(T)-box
	x1=linspace(-5,5,100)
	y1=R*sin(theta)-1./tan(theta)*x1+ R*(cos(theta))**2/sin(theta) -box
	y2=R*sin(theta)+1./tan(theta)*(x1+ R*cos(theta)) -box	
	plot(x1,y1)
	plot(x1,y2)
	xtot=append(x,X)
	ytot=append(y,Y)
	savetxt("input.dat",transpose(array([xtot,ytot])).astype(float32))
	plot(xtot,ytot,lw=2)
	print("il numero di punti  generati is "+str(len(xtot)))
	area(append(xtot,xtot[0]),append(ytot,ytot[0]))



from matplotlib.path import Path
def area(xx,yy):
	'''t=linspace(0,2*pi,1000)
	xx=cos(t)
	yy=sin(t)'''
	myp=Path(transpose([xx,yy]))
	xmax=max(xx)
	ymax=max(yy)
	xmin=min(xx)
	ymin=min(yy)
	size=max([xmax,ymax,abs(xmin),abs(ymin)])
	myx= linspace(-size,size,1000)
	X,Y=meshgrid(myx,myx)
	c=myp.contains_points(transpose(array([X.flatten(),Y.flatten()])))
	imshow(c.reshape(shape(X)))
	dxdy=diff(myx)[0]**2
	print ("area is "+str(sum(c)*dxdy))
