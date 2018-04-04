def velfield(cnum,nop):
	x,y=transpose(loadtxt("posizione/dati_"+str(cnum)))
	vx,vy=transpose(loadtxt("velocity/velocity_"+str(cnum)))
	fx,fy=transpose(loadtxt("forza/forza_"+str(cnum)))
	L=3.
	bns = linspace(-L/2.,L/2.,nop)
	h,bx,by = histogram2d(x,y,bins=bns)
	figure()
	imshow(h,extent=(-L/2,L/2,-L/2,L/2))
	hvx,bx,by = histogram2d(x,y,bins=bns,weights=vx)
	hvy,bx,by = histogram2d(x,y,bins=bns,weights=vy)
	bxx = (bx[1:]+bx[:-1])/2.
	byy = (by[1:]+by[:-1])/2.
	quiver(bxx,byy,hvx,hvy)
	figure()
	imshow(h,extent=(-L/2,L/2,-L/2,L/2))
	hfx,bx,by = histogram2d(x,y,bins=bns,weights=fx)
	hfy,bx,by = histogram2d(x,y,bins=bns,weights=fy)
	bxx = (bx[1:]+bx[:-1])/2.
	byy = (by[1:]+by[:-1])/2.
	quiver(byy,bxx,hfy,hfx,scale=10000)

def velfieldav(imin,imax,nop):
	x,y=transpose(loadtxt("posizione/dati_"+str(imin)))
	vx,vy=transpose(loadtxt("velocity/velocity_"+str(imin)))
	"""fx,fy=transpose(loadtxt("forza/forza_"+str(imin)))
	cnd = (fx!=0)|(fy!=0)
	x = x[cnd]
	y = y[cnd]
	vx = vx[cnd]
	vy = vy[cnd]"""
	L=3.
	bns = linspace(-L/2.,L/2.,nop)
	ha,bx,by = histogram2d(x,y,bins=bns)
	hvxa,bx,by = histogram2d(x,y,bins=bns,weights=vx)
	hvya,bx,by = histogram2d(x,y,bins=bns,weights=vy)
	for i in arange(imin+1,imax+1):	
		print("i = "+str(i)+" / "+str(imax)+"\n")
		x,y=transpose(loadtxt("posizione/dati_"+str(i)))
		vx,vy=transpose(loadtxt("velocity/velocity_"+str(i)))
		"""fx,fy=transpose(loadtxt("forza/forza_"+str(imin)))
		cnd = (fx!=0)|(fy!=0)
		x = x[cnd]
		y = y[cnd]
		vx = vx[cnd]
		vy = vy[cnd]"""
		#L=3.
		#bns = linspace(-L/2.,L/2.,nop)
		h,bx,by = histogram2d(x,y,bins=bns)
		hvx,bx,by = histogram2d(x,y,bins=bns,weights=vx)
		hvy,bx,by = histogram2d(x,y,bins=bns,weights=vy)
		ha += h
		hvxa += hvx
		hvya += hvy
	figure()
	imshow(ha,extent=(-L/2,L/2,-L/2,L/2))
	bxx = (bx[1:]+bx[:-1])/2.
	byy = (by[1:]+by[:-1])/2.
	#hvxa = 0.3*hvxa/abs(hvxa).max()
	#hvya = 0.3*hvya/abs(hvxa).max()
	quiver(bxx,byy,hvya/float(imax-imin),hvxa/float(imax-imin),scale=1000)
	return bxx,byy,hvxa/float(imax-imin),hvya/float(imax-imin)

	
