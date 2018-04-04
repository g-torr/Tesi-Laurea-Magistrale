from pylab import*
import scipy.special
def     xy(r,phi):
        return r*np.cos(phi), r*np.sin(phi)

def	grafica(i=8,b="input.dat"):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	x1,y1=transpose(loadtxt(b))
	x1=append(x1,x1[0])
	y1=append(y1,y1[0])
	plot(x,y,"b.")
	plot(x1,y1,"-g")
	#xlim(-2,2)	
	
def	grafica2(i=8,r1=1.5,r2=0.5):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	N=len(x)
	#plot(x,y,"bo")
	xmn = -r1-0.1
	xmx = r1 +0.1
	ymn = -r2 -0.1
	ymax= r2 +0.1
	h,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,0.05),arange(ymn,ymax,0.05)])	
	h=h/N
	phis=np.arange(0,6.28,0.01)
#	plot( *xy(r,phis), c='r',ls='-' )
	return h,bx,by

def	grafica_all(nof=8,r1=1.5,r2=0.5):
	size=8.
	h,bx,by=grafica2(0,r1,r2)
	for i in arange (1,nof):
		h,bx_par,by_par=grafica2(i,r1,r2)
		bx+=bx_par
		by+=by_par
	h=transpose(h)	
	h=h/nof
	imshow(h,extent=(-r1-0.1,r1+0.1,-r2-0.1,r2+0.1))
	
	

	
def	forza(i,fmax,nb):
	forza= "./forza/forza_"+str(i)
	fx,fy=transpose(loadtxt(forza))
	N=len(fx)
	modulo=sqrt(fx**2+fy**2)
	fx=fx[fx!=0]
	fy=fy[fy!=0]
	part_on_bound=len(modulo[modulo!=0])
	'''print(" le pallette che sbattono sono "+str(part_on_bound)+" su "+str(N))
	print(" forza risultante F_x="+str(sum(fx)))
	print(" forza risultante F_y="+str(sum(fy)))'''
	figure(3)
	plot(i,mean(modulo),"o")
	h,b = histogram(modulo[modulo!=0], bins=linspace(-0.1,fmax,nb),normed=True)
	'''hx,b1 = histogram(fx, bins=linspace(-fmax,fmax,nb),normed=True)
	hy,b1 = histogram(fy, bins=linspace(-fmax,fmax,nb),normed=True)'''
	figure(0)
	plot(b[:-1],h)
	'''figure(1)
	plot(b1[:-1],hx)
	figure(2)
	plot(b1[:-1],hy)'''
	ratio= 1.*part_on_bound/(1.*N)
	#return h,hx,hy,ratio
	return h,ratio
	
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

def	media_temporale(start,finish,r1=1,r2=1,D=1.,tau=1.,nb=150):
	figure(0)
	title("istogramma del modulo  della forza diversa da 0")
	fmax=5*sqrt(D/tau)
	#h_sum,hx_sum,hy_sum,parti_on_bound = forza(0,fmax,nb)
	h_sum,parti_on_bound=forza(start,fmax,nb)
	p_on_b=[]
	p_on_b=append(p_on_b,parti_on_bound)
	for i in arange(start+1,finish+1):
		h,parti_on_bound=forza(i,fmax,nb)
		p_on_b=append(p_on_b,parti_on_bound)
		h_sum+=h
	'''	hx_sum+=hx
		hy_sum+=hy'''
	print("simulated ratio particles on boundary over all particles"+str(mean(p_on_b)))
	theoretical=(normalization(r1,r2,D,tau)-pi*r1*r2)/normalization(r1,r2,D,tau)
	print("theoretical probabilty  particles on boundary"+str(theoretical)) 
	#h_sum=array(h_sum)
	figure(3)
	xlabel("t")
	ylabel("$<|F|>$")
	savefig("modulo_thermalization.pdf")
	modulo_forza(start,finish,r1,r2,D,tau)
	b=linspace(-0.1,fmax,nb)
	f = (b[1:]+b[:-1])/2. #credo serva per prendere il centro
	figure(0)
	plot(f,h_sum,'-k',lw=2)
	b1=linspace(-fmax,fmax,nb)
	'''fx = (b1[1:]+b1[:-1])/2.
	figure(1)
	title("histogram of the x axis component of the force")
	plot(fx,hx_sum,'-k',lw=2)
	figure(2)
	plot(fx,hy_sum,'-k',lw=2)'''	

def	campo_vettoriale(i=8,r1=1.5,r2=0.5):
	x,y=transpose(loadtxt("./posizione/dati_"+str(i)))
	fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
	#modulo=sqrt(fx*fx+fy*fy)
	xmn = -r1-0.1
        xmx = r1 +0.1
        ymn = -r2 -0.1
        ymax= r2 +0.1
	xstep=(xmx-xmn)/100
	ystep=(ymax-ymn)/100
	h,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)])
	hfx,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=fx)
	hfy,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=fy)
	hfx = where(h>0,hfx/h,0)
	hfy = where(h>0,hfy/h,0)
	f=hfx**2+hfy**2
	print(len(bx),len(by))
	quiver(bx[:-1]+xstep,by[:-1]+ystep,hfy,hfx,f,scale=20)
#        quiver(by[:-1],bx[:-1],hfy,hfx,f,scale=20)

	#print("la media del modulo della forza e' "+str(mean(modulo)))
	ax=gca()#fissa gli assi
	ax.set_aspect(1)#mette la stessa scala  tra x e y
	show()
	colorbar()	

def	modulo_forza(start,finish,r1,r2,D,tau):
	fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))
	n=len(fx[(fx!=0)|(fy!=0)])
	somma_x=sum(fx)
	somma_y=sum(fy)
	N=len(fx)
	modulo =sum(sqrt(fx**2+fy**2))
        for i in arange(start+1,finish+1):
		fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_mod=sqrt(fx**2+fy**2)
		modulo+=sum(temp_mod)
		somma_x+=sum(fx)
		somma_y+=sum(fy)
		n+=len(fx[(fx!=0)|(fy!=0)])
		N+=len(fx)
	print("media rispetto alle sole particelle al bordo")
	print("fx medio"+str(somma_x/float(n)))
	print("fy medio"+str(somma_y/float(n)))
	print("media rispetto a tutte")
	print("fx medio"+str(somma_x/float(N)))
	print("fy medio"+str(somma_y/float(N)))	
	print("simulated average of the force module (rispetto a tutte)"+str(modulo/float(N)))	
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	theoretical=(D*arc_lenght+ D*sqrt(D*tau)*2.*pi	)/normalization(r1,r2,D,tau)
	print(" theoretical average force module( rispetto a tutte)"+str(theoretical))
	print("simulated conditional  average of the force module (to the boundary)"+str(modulo/float(n)))


def f(x,A,B):
	return A*x + B
from scipy.optimize import curve_fit

def	probability_in_curvature(r1,r2,D,tau,nof=9):
	x=[]
	y=[]
	fx=[]
	fy=[]
	for i in arange(0,nof):
		x_temp,y_temp=transpose(loadtxt("./posizione/dati_"+str(i)))
		x=append(x,x_temp)
		y=append(y,y_temp)
		fx_temp,fy_temp=transpose(loadtxt("./forza/forza_"+str(i)))
		fx=append(fx,fx_temp)
		fy=append(fy,fy_temp)
	'''xmn=min(x)
	xmx=max(x)
	ymn=min(y)
	ymax=max(y)'''
	app_x=[]
	app_y=[]
	modulo=sqrt(fx**2+fy**2)
	for i in arange(0,len(x)):
		if(modulo[i]!=0):
			app_x=append(app_x,x[i])
			app_y=append(app_y,y[i])
	x,y=transpose(loadtxt("input.dat"))
	count=[]
	x=append(x,x[0])
	y=append(y,y[0])
	for i in arange(0,(len(x)-1)):
		c=0
		for j in arange(0,len(app_x)):
			if ((min(x[i],x[i+1])<app_x[j]<max(x[i],x[i+1]))&(min(y[i],y[i+1])<app_y[j]<max(y[i],y[i+1]))):
		 		c+=1
		count=append(count,c/sqrt((x[i]-x[i+1])**2+(y[i]-y[i+1])**2))
		#count=append(count,c)
	k=np.load("curvatura.npy")	
	count0=count/len(modulo)
	count1=count/len(modulo[modulo!=0])
	plot(k[:-1],count,".")
	np.save("probability_in_curv.npy",array([k[:-1],count0]))
	np.save("probability_in_curv_cond.npy",array([k[:-1],count1]))
	A,B= curve_fit(f, k[:-1],count)[0]
	plot(k[:-1],f(k[:-1],A,B),"-")	
	print("interpolation with y= "+str(A)+" x + "+str(B))
	print("theoretical line y="+str(D*tau/normalization(r1,r2,D,tau))+" x+"+str(sqrt(D*tau)/normalization(r1,r2,D,tau)))
	return A,B

	'''	
	h,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)])
	hf,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=modulo)
	'''	
def normalization(r1,r2,D,tau):
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	return 	pi*r1*r2+sqrt(D*tau)*arc_lenght+D*tau*2*pi	
	
		
