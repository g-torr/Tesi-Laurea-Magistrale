def	grafica(i=8,b="input2.dat"):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	x1,y1=transpose(loadtxt(b))
	x1=append(x1,x1[0])
	y1=append(y1,y1[0])
	plot(x,y,"b.")
	plot(x1,y1,"-g")
	#xlim(-2,2)	
	
def	grafica2(i=8,r=1.,size=3.):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	N=len(x)
	xmn = -size/2.
	xmx = size/2.
	h,bx,by = histogram2d(x,y,bins=arange(xmn,xmx,0.03))	
	h=h/N
	phis=np.arange(0,6.28,0.01)
#	plot( *xy(r,phis), c='r',ls='-' )
	return h,bx,by

def	grafica_all(nof=8):
	size=3.
	r=1.
	h,bx,by=grafica2(i=0)
	for i in arange (1,nof):
		h,bx_par,by_par=grafica2(i,r,size)
		bx+=bx_par
		by+=by_par
	h=transpose(h)	
	h=h/nof
	imshow(h,extent=(-size/2.,size/2.,-size/2.,size/2.))
	
	
	
def	forza(i,fmax,nb):
	forza= "./forza/forza_"+str(i)
	fx,fy=transpose(loadtxt(forza))
	N=len(fx)
	fx=fx[fx!=0]
	print(" le pallette che sbattono sono "+str(len(fx))+" su "+str(N))#e' vero se il poligono scelto non ha tratti ne' orizzontali ne' verticali
	print(" forza risultante F_x="+str(sum(fx)))
	print(" forza risultante F_y="+str(sum(fy)))
	
	h,b = histogram(fx, bins=linspace(-fmax,fmax,nb),normed=True)
	plot(b[:-1],h)	
	return h	
	#hist(fx,10,normed=True)
	
def	momenti(i=8):
	forza= "./forza/forza_"+str(i)
	fx,fy=transpose(loadtxt(forza))
	N=len(fx)
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	mom=[]
	for i in range(0,N):
		momento=x[i]*fy[i]-y[i]*fx[i]
		mom+=[momento]	
	mom=array(mom)	
	print("il momento risultante rispetto all'origine="+str(sum(mom)))	
	mom1=mom[mom!=0]
	print("i momenti diveri da 0 sono"+str(len(mom1)))
	# ***************plot i momenti**********
	title("istogramma dei momenti diversi da 0")
	hist(mom1,20,normed=True)
	'''h,b = histogram(mom1, bins=linspace(mom1.min(),mom1.max(),100),normed=True)
	x = (b[1:]+b[:-1])/2. #credo serva per predendere il centro
	plot(x,h,'-ob')
	'''

def	media_temporale(nb):
	title("istogramma della F_x diversa da 0")
	D=1	
	tau=0.06
	fmax=4*sqrt(D/tau)
	h_sum = forza(0,fmax,nb)
	for i in arange(1,8):
		h=forza(i,fmax,nb)
		h_sum+=h
	#h_sum=array(h_sum)
	b=linspace(-fmax,fmax,nb)
	x = (b[1:]+b[:-1])/2. #credo serva per predendere il centro
	plot(x,h_sum,'-k',lw=2)
	modulo_forza()
'''	print("len(h)="+str(len(h)))
	print("len(b)="+str(len(b)))'''
def	campo_vettoriale(i=8,size=3.):
	x,y=transpose(loadtxt("./posizione/dati_"+str(i)))
	fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
	h,bx,by = histogram2d(x,y,bins=arange(-size/2.,size/2,.1))
	hfx,bx,by = histogram2d(x,y,bins=arange(-size/2.,size/2.,.1),weights=fx)
	hfy,bx,by = histogram2d(x,y,bins=arange(-size/2,size/2,.1),weights=fy)
	hfx = where(h>0,hfx/h,0)
	hfy = where(h>0,hfy/h,0)
	f=sqrt(hfx**2+hfy**2)
	quiver(by[:-1],bx[:-1],hfx,hfy,f,scale=20)
#        quiver(by[:-1],bx[:-1],hfy,hfx,f,scale=20)
	ax=gca()#fissa gli assi
	ax.set_aspect(1)#mette la stessa scala  tra x e y
	show()
	colorbar()	

def	modulo_forza():
	D=1.
	tau=0.06
	fx,fy=transpose(loadtxt("./forza/forza_0"))
	modulo =sqrt(fx**2+fy**2)
	for i in arange(1,8):
		fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_mod=sqrt(fx**2+fy**2)
		append(modulo,temp_mod)
	print("la media del modulo della forza e' "+str(mean(modulo)))
	print("La UCNA prevede che il modulo della forza ="+str(sqrt(D*pi/(2*tau))))
