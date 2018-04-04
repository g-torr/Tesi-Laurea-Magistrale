def 	xy(r,phi):
  	return r*np.cos(phi), r*np.sin(phi)

def	grafica(i=8,r=1.):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	plot(x,y,"bo") 
	phis=np.arange(0,6.28,0.01)
	plot( *xy(r,phis), c='r',ls='-' )	
	#xlim(-2,2)	
	
def	grafica2(i=8,r=1.):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	xmn = x.min()
	xmx = x.max()
	h,bx,by = histogram2d(x,y,bins=arange(xmn,xmx,0.05))	
	h=transpose(h)	
	imshow(h,extent=(xmn,xmx,xmn,xmx))
	#ora plotto la circonferenza	
	phis=np.arange(0,6.28,0.01)
	plot( *xy(r,phis), c='r',ls='-' )
	#plot(x1,y1,"-go")


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

def	modulo_forza(i):
	forza= "./forza/forza_"+str(i)
	fx,fy=transpose(loadtxt(forza))
#	fx=fx[fx!=0]	
#	fy=fy[fy!=0]
	modulo=fx*fx+fy+fy
	return mean(modulo)

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

def	media_temporale(nb):#nb number of bins
	title("istogramma della F_x diversa da 0")
	D=0.1	
	tau=0.6
	fmax=4*sqrt(D/tau)
	h_sum = forza(0,fmax,nb)
	f_sum=0	
	for i in arange(1,8):
		h=forza(i,fmax,nb)
		h_sum+=h
		f_sum+=modulo_forza(i)
	b=linspace(-fmax,fmax,nb)
	x = (b[1:]+b[:-1])/2. #credo serva per predendere il centro e perche linspace mi da un array con un elemento in piu
	plot(x,h_sum,'-k',lw=2)#plot in black, line width=2
	print("la somma dei quadrati dei moduli ="+str(f_sum))
