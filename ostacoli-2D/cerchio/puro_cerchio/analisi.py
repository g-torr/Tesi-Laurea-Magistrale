def     xy(r,phi):
        return r*np.cos(phi), r*np.sin(phi)

def	grafica(i=8,b="input2.dat"):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	x1,y1=transpose(loadtxt(b))
	x1=append(x1,x1[0])
	y1=append(y1,y1[0])
	plot(x,y,"b.")
	plot(x1,y1,"-g")
	#xlim(-2,2)	
	
def	grafica2(i=8,size=3.,r=1.):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	#plot(x,y,"bo")
	N=len(x)
	xmn = -size/2.
	xmx = size/2.
	h,bx,by = histogram2d(x,y,bins=arange(xmn,xmx,0.05))	
	h=h/N
	phis=np.arange(0,6.28,0.01)
	plot( *xy(r,phis), c='r',ls='-' )
	return h,bx,by

def	grafica_all(nof=8,r=1.):
	size=3.
	h,bx,by=grafica2(i=0)
	for i in arange (1,nof):
		h,bx_par,by_par=grafica2(i,size,r)
		bx+=bx_par
		by+=by_par
	h=transpose(h)	
	h=h/nof
	imshow(h,extent=(-size/2.,size/2.,-size/2.,size/2.))
	
	
	
def	forza(i,fmax,nb):
	forza= "./forza/forza_"+str(i)
	fx,fy=transpose(loadtxt(forza))
	N=len(fx)
	modulo=sqrt(fx**2+fy**2)
	modulo=modulo[modulo!=0]
	fx=fx[fx!=0]
	fy=fy[fy!=0]
	part_on_bound=len(modulo)
	print(" le pallette che sbattono sono "+str(part_on_bound)+" su "+str(N))
	print(" forza risultante F_x="+str(sum(fx)))
	print(" forza risultante F_y="+str(sum(fy)))
	h,b = histogram(modulo, bins=linspace(-0.1,fmax,nb),normed=True)
	hx,b1 = histogram(fx, bins=linspace(-fmax,fmax,nb),normed=True)
	hy,b1 = histogram(fy, bins=linspace(-fmax,fmax,nb),normed=True)
	figure(0)
	plot(b[:-1],h)
	figure(1)
	plot(b1[:-1],hx)
	figure(2)
	plot(b1[:-1],hy)
	ratio= 1.*part_on_bound/(1.*N)
	return h,hx,hy,ratio
	
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

def	media_temporale(nb=150):
	figure(0)
	title("istogramma del modulo  della forza diversa da 0")
	D=0.1	
	tau=0.06
	fmax=5*sqrt(D/tau)
	h_sum,hx_sum,hy_sum,parti_on_bound = forza(0,fmax,nb)
	p_on_b=[]
	p_on_b=append(p_on_b,parti_on_bound)
	for i in arange(1,8):
		h,hx,hy,parti_on_bound=forza(i,fmax,nb)
		p_on_b=append(p_on_b,parti_on_bound)
		h_sum+=h
		hx_sum+=hx
		hy_sum+=hy
	print("simulated ratio particles on boundary over all particles"+str(mean(p_on_b)))
	theoretical=2.*pi*sqrt(D*tau)*(1.+sqrt(D*tau))/(pi+2*pi*sqrt(D*tau))*(1+sqrt(D*tau))
	print("theoretical probabilty  particles on boundary"+str(theoretical)) 
	#h_sum=array(h_sum)
	modulo_forza()
	print("theoretical one"+str(2*D*(1+sqrt(D*tau))/(1+2*D*tau+2*sqrt(D*tau))))
	b=linspace(-0.1,fmax,nb)
	f = (b[1:]+b[:-1])/2. #credo serva per predendere il centro
	figure(0)
	plot(f,h_sum,'-k',lw=2)
	b1=linspace(-fmax,fmax,nb)
	fx = (b1[1:]+b1[:-1])/2.
	figure(1)
	title("histogram of the x axis component of the force")
	plot(fx,hx_sum,'-k',lw=2)
	figure(2)
	plot(fx,hy_sum,'-k',lw=2)		

'''	print("len(h)="+str(len(h)))
	print("len(b)="+str(len(b)))'''
def	campo_vettoriale(i=8,size=3.):
	x,y=transpose(loadtxt("./posizione/dati_"+str(i)))
	fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
	#modulo=sqrt(fx*fx+fy*fy)
	h,bx,by = histogram2d(x,y,bins=arange(-size/2.,size/2,.1))
	hfx,bx,by = histogram2d(x,y,bins=arange(-size/2.,size/2.,.1),weights=fx)
	hfy,bx,by = histogram2d(x,y,bins=arange(-size/2,size/2,.1),weights=fy)
	hfx = where(h>0,hfx/h,0)
	hfy = where(h>0,hfy/h,0)
	f=hfx**2+hfy**2
	quiver(by[:-1]+0.05,bx[:-1],hfy,hfx,f,scale=20)
#        quiver(by[:-1],bx[:-1],hfy,hfx,f,scale=20)

	#print("la media del modulo della forza e' "+str(mean(modulo)))
	ax=gca()#fissa gli assi
	ax.set_aspect(1)#mette la stessa scala  tra x e y
	show()
	colorbar()	
def	modulo_forza():
	fx,fy=transpose(loadtxt("./forza/forza_0"))
	modulo =sqrt(fx**2+fy**2)
	for i in arange(1,8):
		fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_mod=sqrt(fx**2+fy**2)
		append(modulo,temp_mod)
	print("la media del modulo della forza e' "+str(mean(modulo)))

