def mkframe(cnum,nob):
	L = 3.
	bns = linspace(-L/2,L/2.,nob)
	x,y = transpose(loadtxt("posizione/dati_"+str(cnum)))
	vx,vy = transpose(loadtxt("velocity/velocity_"+str(cnum)))
	"""x = x[0:nop]
	y = y[0:nop]
	vx = vx[0:nop]
	vy = vy[0:nop]"""
	"""print vx
	print vy"""
	h,bx,by = histogram2d(x,y,bins=bns)
	"""figure()
	imshow(transpose(h),extent=(-L/2.,L/2.,-L/2.,L/2.),origin="lower")	
	quiver(x,y,vx,vy)
	plot(x,y,"ok")"""
	#figure()
	#imshow(transpose(h),extent=(-L/2.,L/2.,-L/2.,L/2.),origin="lower")	
	hvx,bx,by = histogram2d(x,y,bins=bns,weights=vx)
	hvy,bx,by = histogram2d(x,y,bins=bns,weights=vy)
	hx,bx,by = histogram2d(x,y,bins=bns,weights=x)
	hy,bx,by = histogram2d(x,y,bins=bns,weights=y)
	#h,bx,by = histogram2d(x,y,bins=bns)
	h = h.flatten()
	hx = hx.flatten()
	hy = hy.flatten()
	hvx = hvx.flatten()
	hvy = hvy.flatten()
	hx = compress(h>0,hx)/compress(h>0,h)
	hy = compress(h>0,hy)/compress(h>0,h)
	hJx = compress(h>0,hvx)
	hJy = compress(h>0,hvy)
	hvx=hJx/compress(h>0,h)
	hvy=hJy/compress(h>0,h)
	#quiver(hx,hy,hvx,hvy,color="w")
	"""print hvx[hvx!=0]
	print hvy[hvy!=0]"""
	return h,hx,hy,hvx,hvy,hJx,hJy	

def media(start,finish,nob=100):
	h_sum,hx,hy,hvx_sum,hvy_sum,hJx_sum,hJy_sum= mkframe(start,nob)
	for i in arange(start+1,finish+1):
		print ("siamo a "+str(i)+" su "+str(finish)+"\n")
		h,hx,hy,hvx,hvy,hJx,hJy= mkframe(i,nob)
		h_sum	+=h
		hvx_sum +=hvx
		hvy_sum +=hvy
		hJx_sum +=hJx
		hJy_sum +=hJy
	print("sum(Jx)= "+str(sum(hJx_sum)/float(sum(h_sum)))+"sum(Jy)= "+str(sum(hJy_sum)/float(sum(h_sum))))	
	print("sum(Vx)= "+str(sum(hvx_sum)/float(sum(h_sum)))+"sum(Jy)= "+str(sum(hvy_sum)/float(sum(h_sum))))
	quiver(hx,hy,hvx_sum/float(finish-start),hvy_sum/float(finish-start),color="w")
	return hx,hy,hvx_sum/float(finish-start),hvy_sum/float(finish-start),hJx_sum/float(finish-start),hJy_sum/float(finish-start)
	

		
		 
def trsp(x):
	return transpose(x)
