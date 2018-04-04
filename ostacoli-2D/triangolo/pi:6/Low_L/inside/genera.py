from matplotlib.path import Path

def genera(theta=pi/6):
	x=[cos(theta)	,	0,		-cos(theta),	0]
	y=[sin(theta)-1.,1/sin(theta)-1,sin(theta)-1.,		0]
	plot(x,y)
	figure()
	savetxt("input.dat",transpose(array([x,y])).astype(float32))
	print("il numero di vertici = "+str(len(x)))
	surf=area(x,y)
	np.save("area.npy",[surf,3])
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
	return sum(c)*dxdy

