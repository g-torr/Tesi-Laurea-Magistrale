def all_analisi():
	tau=[0.3,0.5,0.06,1,2,3,4,10,20]
	J=[]
	for t in tau:
		fonam="tau="+str(t)
		J+=[media_modulo_corrente(fonam)]
		print "tau="+str(t)
	plot(tau,J,"-o")
	np.save("media_J.npy",array([tau,J]))

def mkframe(fonam,cnum,nob):
	L = 3.
	bns = linspace(-L/2,L/2.,nob)
	x,y = transpose(loadtxt(fonam+"/posizione/dati_"+str(cnum)))
	vx,vy = transpose(loadtxt(fonam+"/velocity/velocity_"+str(cnum)))
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
	return hx,hy,hvx,hvy,hJx,hJy	

def media(fonam,start,finish,nob=100):
	hx,hy,hvx_sum,hvy_sum,hJx_sum,hJy_sum= mkframe(fonam,start,nob)
	for i in arange(start+1,finish+1):
		print ("siamo a "+str(i)+" su "+str(finish)+"\n")
		hx,hy,hvx,hvy,hJx,hJy= mkframe(fonam,i,nob)
		hvx_sum +=hvx
		hvy_sum +=hvy
		hJx_sum +=hJx
		hJy_sum +=hJy
	quiver(hx,hy,hvx_sum/float(finish-start),hvy_sum/float(finish-start),color="w")
	return hx,hy,hvx_sum/float(finish-start),hvy_sum/float(finish-start),hJx_sum/float(finish-start),hJy_sum/float(finish-start)
	

		
		 
def media_modulo_corrente(fonam):
	hx,hy,hvx,hvy,hJx,hJy=media(fonam,18,19)
	hJ=sqrt(hJx**2+hJy**2)	#media del modulo della corrente
	return mean(hJ)
