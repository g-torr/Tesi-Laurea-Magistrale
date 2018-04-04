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
	
def	grafica2(i,xmin,xmax,ymin,ymax):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	N=len(x)
	#plot(x,y,"bo")
	h,bx,by = histogram2d(x,y,bins=[arange(xmin,xmax,0.02),arange(ymin,ymax,0.02)])	
	h=h/N
	phis=np.arange(0,6.28,0.01)
#	plot( *xy(r,phis), c='r',ls='-' )
	return h,bx,by

def	grafica_all(nof=8,r1=1.5,r2=0.5):
	size=8.
	posizione= "./posizione/dati_0"
	x,y=transpose(loadtxt(posizione))
	xmin=min(x)-0.1
	ymin=min(y)-0.1
	xmax=max(x)+0.1
	ymax=max(y)+0.1
	h,bx,by=grafica2(0,xmin,xmax,ymin,ymax)
	for i in arange (1,nof):
		h,bx_par,by_par=grafica2(i,xmin,xmax,ymin,ymax)
		bx+=bx_par
		by+=by_par
	h=transpose(h)	
	h=h/nof
	imshow(h,extent=(xmin,xmax,ymin,ymax))
	
def	campo_velocity(nf=8):
	x,y=transpose(loadtxt("./posizione/dati_"+str(0)))
	vx,vy=transpose(loadtxt("./velocity/velocity_"+str(0)))
	somma_vx=sum(vx)
	somma_vy=sum(vy)
	N=len(vx)
	for i in arange (1,nf):					
		x,y=transpose(loadtxt("./posizione/dati_"+str(i)))
		vx,vy=transpose(loadtxt("./velocity/velocity_"+str(i)))
		somma_vx+=sum(vx)
		somma_vy+=sum(vy)
		N+=len(vx)
	print("somma delle vy"+str(somma_vx/float(N)))
	print("somma delle vx"+str(somma_vy/float(N)))
	

def	forza(i,fmax,nb):
	forza= "./forza/forza_"+str(i)
	fx,fy=transpose(loadtxt(forza))
	N=len(fx)
	modulo=sqrt(fx**2+fy**2)
	modulo=modulo[modulo!=0]
	fx=fx
	fy=fy
	part_on_bound=len(modulo)
	#print(" le pallette che sbattono sono "+str(part_on_bound)+" su "+str(N))
	'''print(" forza media F_x="+str(mean(fx)))
	print(" forza media F_y="+str(mean(fy)))'''
	figure(4)
	plot(i,mean(fy),"o")
	figure(3)
	plot(i,mean(fx),"o")
        hx,bx = histogram(fx[fx!=0], bins=linspace(-fmax,fmax,nb),normed=True)
	hy,by = histogram(fy[fy!=0], bins=linspace(-fmax,fmax,nb),normed=True)
	'''hx,b1 = histogram(fx, bins=linspace(-fmax,fmax,nb),normed=True)
	hy,b1 = histogram(fy, bins=linspace(-fmax,fmax,nb),normed=True)'''
	bx = (bx[1:]+bx[:-1])/2. #credo serva per prendere il centro
        by = (by[1:]+by[:-1])/2.
	figure(0)
	title("fx")
	plot(bx,hx)
	figure(1)
	title("fy")
	plot(by,hy)
	ratio= float(part_on_bound)/float(N)
	#return h,hx,hy,ratio
	return hx,hy,bx,by,ratio
	
def	momenti(i=8,bins=50):
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
	print("i momenti diversi da 0 sono"+str(len(mom1)))
	# ***************plot i momenti**********
	title("istogramma dei momenti diversi da 0")
	hist(mom1,bins,normed=True)
	'''h,b = histogram(mom1, bins=linspace(mom1.min(),mom1.max(),100),normed=True)
	x = (b[1:]+b[:-1])/2. #credo serva per predendere il centro
	plot(x,h,'-ob')
	'''

def	media_temporale(start,stop,nb=150):
	figure(0)
	title("istogramma del modulo  della forza diversa da 0")
	D=50.	
	tau=10.
	fmax=3*sqrt(D/tau)
	#h_sum,hx_sum,hy_sum,parti_on_bound = forza(0,fmax,nb)
	hx,hy,bx,by,parti_on_bound=forza(start,fmax,nb)
	p_on_b=[]
	p_on_b=append(p_on_b,parti_on_bound)
	for i in arange(start+1,stop+1):
		par_x,par_y,bx,by,parti_on_bound=forza(i,fmax,nb)
		p_on_b=append(p_on_b,parti_on_bound)
		hx+=par_x
		hy+=par_y
	print("simulated ratio particles on boundary over all particles "+str(mean(p_on_b)))
	#theoretical=(normalization(r1,r2)-pi*r1*r2)/normalization(r1,r2)
	#print("theoretical probabilty  particles on boundary"+str(theoretical)) 
	#h_sum=array(h_sum)
	modulo_forza(start,stop)
	'''b=linspace(-0.1,fmax,nb)
	bx = (bx[1:]+bx[:-1])/2. #credo serva per prendere il centro
	by = (by[1:]+by[:-1])/2.
	figure(0)
	xlabel("fx")
	ylabel("density")
	plot(bx,hx,'-k',lw=2)'''
	figure(1)
        xlabel("fy")
        ylabel("density")
	figure(4)
	title ("$F_y$ component")
	savefig("Fy_termalization.pdf")
	figure(3)
	savefig("Fx_termalization.pdf")
	title ("$F_x$ component")
	figure()
	plot(by,hy,'-k',lw=2)
	#b1=linspace(-fmax,fmax,nb)
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
def forza2(start,stop):
	area,arco=transpose(loadtxt("area.txt"))
        fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))
	x,y=transpose(loadtxt("./posizione/dati_"+str(start)))
	N=len(fy)
	fy=fy
	for i in arange(start+1,stop+1):
                temp_fx,temp_fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_x,temp_y=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx=append(fx,temp_fx)
		x=append(x,temp_x)
		fy=append(fy,temp_fy)
		y=append(y,temp_y)
	r0=0.7
	r=sqrt(x**2+y**2)
	rho=count_nonzero(r<r0)/(pi*r0**2) #densita' di bulk
	h,b=histogram(y,bins=arange(-1.,max(y)+0.02,0.0005),weights=fy) #it is the yforce density distribution, density is normalized at N
	h=h/rho/(stop-start+1)#normalizzo rispetto alla densita' di bulk
	sh=cumsum(h)
	b=(b[:-1]+b[1:])/2
	figure(0)
	title ("$F_y(y)$")
	xlabel("$y$")
	ylabel ("$F_y$")
	np.save("F.npy",array([b,h]))
	plot(b,h)
	figure(1)
	title ("cumulative $F_y(y)$")
	xlabel("$y$")
	ylabel ("$cum F_y$")
	plot(b,sh)
	np.save("F_cumulative.npy",array([b,sh]))
			
def	forza_netta(nf=9):
        fx,fy=transpose(loadtxt("./forza/forza_0"))
        nx=len(fx[fx!=0])
	ny=len(fy[fy!=0])
	somma_x=sum(fx)
	somma_y=sum(fy)
	N=len(fx)
        for i in arange(1,nf):
                tempx,tempy=transpose(loadtxt("./forza/forza_"+str(i)))
		somma_x+=sum(tempx)
		somma_y+=sum(tempy)
		nx+=len(tempx[tempx!=0])
		ny+=len(tempy[tempy!=0])
		N+=len(tempx)
	print("viene calcolata mediando solo sulle particelle al bordo")
	print("fx medio"+str(somma_x/float(nx)))
	print("fy medio"+str(somma_y/float(ny)))
	print("viene calcolata mediando su tutte")
	print("fx medio"+str(somma_x/float(N)))
	print("fy medio"+str(somma_y/float(N)))

def	modulo_forza(start,stop):
	fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))
	figure()
	modulo =sqrt(fx**2+fy**2)
	plot(start,mean(modulo),"o")
	for i in arange(start+1,stop+1):
		fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_mod=sqrt(fx**2+fy**2)
		plot(i,mean(temp_mod),"o")
		append(modulo,temp_mod)
	print("la media del modulo della forza e' "+str(mean(modulo)))
	#print(" theoretical average force module"+str(theoretical))
	ylabel("module Force")
	savefig("module_thermalization.pdf")
def f(x,A,B):
	return A*x + B
from scipy.optimize import curve_fit

def	probability_in_curvature(r1,r2,nof=9):
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
	D=0.1
	tau=0.06 	
	print("interpolation with y= "+str(A)+" x + "+str(B))
	print("theoretical line y="+str(D*tau/normalization(r1,r2))+" x+"+str(sqrt(D*tau)/normalization(r1,r2)))
	return A,B

	'''	
	h,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)])
	hf,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=modulo)
	'''	
def normalization(r1,r2):
	D=0.1
	tau=0.06
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	return 	pi*r1*r2+sqrt(D*tau)*arc_lenght+D*tau*2*pi	
	
		
