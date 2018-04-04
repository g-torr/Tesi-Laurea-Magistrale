def mkframe(cnum,nob,nop):
	L = 3.
	bns = linspace(-L/2,L/2.,nob)
	x,y = transpose(loadtxt("posizione/dati_"+str(cnum)))
	vx,vy = transpose(loadtxt("velocity/velocity_"+str(cnum)))
	x = x[0:nop]
	y = y[0:nop]
	vx = vx[0:nop]
	vy = vy[0:nop]
	print vx
	print vy
	h,bx,by = histogram2d(x,y,bins=bns)
	imshow(transpose(h),extent=(-L/2.,L/2.,-L/2.,L/2.),origin="lower",interpolation="nearest")	
	quiver(x,y,vx,vy)
	hvx,bx,by = histogram2d(x,y,bins=bns,weights=vx)
	hvy,bx,by = histogram2d(x,y,bins=bns,weights=vy)
	hx,bx,by = histogram2d(x,y,bins=bns,weights=x)
	hy,bx,by = histogram2d(x,y,bins=bns,weights=y)
	h,bx,by = histogram2d(x,y,bins=bns)
	h = h.flatten()
	hx = hx.flatten()
	hy = hy.flatten()
	hvx = hvx.flatten()
	hvy = hvy.flatten()
	hx = compress(h>0,hx)/compress(h>0,h)
	hy = compress(h>0,hy)/compress(h>0,h)
	hvx = compress(h>0,hvx)
	hvy = compress(h>0,hvy)
	quiver(hx,hy,hvx,hvy,color="w",scale=50)
	print hvx[hvx!=0]
	print hvy[hvy!=0]
	plot(x,y,"ok")

def trsp(x):
	return transpose(x)
