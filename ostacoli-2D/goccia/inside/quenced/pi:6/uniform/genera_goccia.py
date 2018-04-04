def genera(R=3,r=0.5,theta=pi/8,n=50):
#	t=linspace(theta,pi-theta,n)      #arco piccolo
	T=linspace((pi-theta),2*pi+theta,n)	#arco grande
#	box=(R-r)/sin(theta)/2.
	box=0
	x=0
	y=R/sin(theta)-box
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
	normalizzazione(R,r,theta)
	area(append(xtot,xtot[0]),append(ytot,ytot[0]))
def normalizzazione(R,r,theta):
	l=(R-r)/tan(theta)
	arc_length=R*(pi+2 *theta)+r*(pi-2*theta)+2*l
	area=2*(r+R)*l/2+R**2*(pi+2*theta)/2+r**2*(pi-2*theta)/2
	savetxt("area.txt",[area,arc_length])

	print("arc_lenght = "+str(arc_length))	
	print("area = "+str(area))		



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
