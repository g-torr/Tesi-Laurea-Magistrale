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
#	print("area is"+str(sum(c)/float(len(myx)**2)))	
