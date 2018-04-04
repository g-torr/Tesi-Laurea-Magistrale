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
	
def	grafica2(i,xmin,xmax,ymin,ymax,pixel):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	N=len(x)
	#plot(x,y,"bo")
	h,bx,by = histogram2d(x,y,bins=[arange(xmin,xmax,pixel),arange(ymin,ymax,pixel)])	
	h=h/N
	phis=np.arange(0,6.28,0.01)
#	plot( *xy(r,phis), c='r',ls='-' )
	return h,bx,by

def	grafica_all(nof=8,pixel=0.02):
	size=8.
	posizione= "./posizione/dati_0"
	x,y=transpose(loadtxt(posizione))
	xmin=min(x)-0.1
	ymin=min(y)-0.1
	xmax=max(x)+0.1
	ymax=max(y)+0.1
	h,bx,by=grafica2(0,xmin,xmax,ymin,ymax,pixel)
	for i in arange (1,nof):
		h,bx_par,by_par=grafica2(i,xmin,xmax,ymin,ymax,pixel)
		bx+=bx_par
		by+=by_par
	h=transpose(h)	
	h=h/nof
	imshow(h,extent=(xmin,xmax,ymin,ymax))
	return h	
	

	
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
	print(" forza media F_x="+str(mean(fx)))
	print(" forza media F_y="+str(mean(fy)))
        hx,bx = histogram(fx, bins=linspace(-fmax,fmax,nb),normed=True)
	hy,by = histogram(fy, bins=linspace(-fmax,fmax,nb),normed=True)
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
	ratio= 1.*part_on_bound/(1.*N)
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

def	media_temporale(nf=9,nb=150):
	figure(0)
	title("istogramma del modulo  della forza diversa da 0")
	D=20.	
	tau=2.
	fmax=4*sqrt(D/tau)
	#h_sum,hx_sum,hy_sum,parti_on_bound = forza(0,fmax,nb)
	hx,hy,bx,by,parti_on_bound=forza(0,fmax,nb)
	p_on_b=[]
	p_on_b=append(p_on_b,parti_on_bound)
	for i in arange(1,nf):
		par_x,par_y,bx,by,parti_on_bound=forza(i,fmax,nb)
		p_on_b=append(p_on_b,parti_on_bound)
		hx+=par_x
		hy+=par_y
	'''	hx_sum+=hx
		hy_sum+=hy'''
	
	print("simulated ratio particles on boundary over all particles"+str(mean(p_on_b)))
	'''theoretical=(normalization(r1,r2)-pi*r1*r2)/normalization(r1,r2)
	print("theoretical probabilty  particles on boundary"+str(theoretical)) 
	#h_sum=array(h_sum)'''
	'''b=linspace(-0.1,fmax,nb)'''
	'''bx = (bx[1:]+bx[:-1])/2. #credo serva per prendere il centro
	by = (by[1:]+by[:-1])/2.'''
	figure(0)
	xlabel("fx")
	ylabel("density")
	plot(bx,hx,'-k',lw=2)
	figure(1)
        xlabel("fy")
        ylabel("density")

	plot(by,hy,'-k',lw=2)
	modulo_forza(D,tau)
	
	'''fx = (b1[1:]+b1[:-1])/2.
	figure(1)
	title("histogram of the x axis component of the force")
	plot(fx,hx_sum,'-k',lw=2)
	figure(2)
	plot(fx,hy_sum,'-k',lw=2)'''	

def	campo_forze(nf):
	x,y=transpose(loadtxt("./posizione/dati_"+str(0)))
	fx,fy=transpose(loadtxt("./forza/forza_"+str(0)))
	#modulo=sqrt(fx*fx+fy*fy)
	xmn = min(x)
        xmx = max(x)
        ymn = min(y)
        ymax= max(y)
	xstep=(xmx-xmn)/100
	ystep=(ymax-ymn)/100
	h,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)])
	hfx,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=fx)
	hfy,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=fy)
	hfx = where(h>0,hfx/h,0)
	hfy = where(h>0,hfy/h,0)
	for i in arange (1,nf):
		h,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)])
		hfx_par,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=fx)
		hfy_par,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=fy)
		hfx_par = where(h>0,hfx_par/h,0)
		hfy_par = where(h>0,hfy_par/h,0)
	hfx+=hfx_par
	hfy+=hfy_par

	f=hfx**2+hfy**2
	print(len(bx),len(by))
	quiver(by[:-1]+ystep,bx[:-1]+xstep,hfy/nf,hfx/nf,f,scale=10)
#        quiver(by[:-1],bx[:-1],hfy,hfx,f,scale=20)

	#print("la media del modulo della forza e' "+str(mean(modulo)))
	ax=gca()#fissa gli assi
	ax.set_aspect(1)#mette la stessa scala  tra x e y
	title("force field")
	show()
	colorbar()	

def	forza_netta(nf=9):
        fx,fy=transpose(loadtxt("./forza/forza_0"))
	Nx=len(fx)
	Ny=len(fy)
        nx=len(fx!=0)
	ny=len(fy!=0)
	somma_x=sum(fx)
	somma_y=sum(fy)
        for i in arange(1,nf):
                tempx,tempy=transpose(loadtxt("./forza/forza_"+str(i)))
		somma_x+=sum(tempx)
		somma_y+=sum(tempy)
		nx+=len(fx[fx!=0])
		ny+=len(fy[fy!=0])
		Nx+=len(fx)
		Ny+=len(fy)	
	'''print("viene calcolata mediando rispetto alle  particelle al bordo")
	print("fx medio"+str(somma_x/nx))
	print("fy medio"+str(somma_y/ny))'''
	print("viene calcolata mediando rispetto a tutte le particelle")
	print("Fx medio"+str(somma_x/Nx))
	print("Fy medio"+str(somma_y/Ny))
def	modulo_forza(D=0.1,tau=0.06):
	fx,fy=transpose(loadtxt("./forza/forza_0"))
	modulo =sqrt(fx**2+fy**2)
	for i in arange(1,8):
		fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_mod=sqrt(fx**2+fy**2)
		append(modulo,temp_mod)
	print("la media del modulo della forza e' "+str(mean(modulo)))
	'''arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	theoretical=(D*arc_lenght+ D*sqrt(D*tau)*2.*pi	)/normalization(r1,r2)
	print(" theoretical average force module"+str(theoretical))'''

def f(x,A,B):
	return A*x + B
from scipy.optimize import curve_fit

def	campo_velocity(nf=8):
	x,y=transpose(loadtxt("./posizione/dati_"+str(0)))
	vx,vy=transpose(loadtxt("./velocity/velocity_"+str(0)))
	#modulo=sqrt(fx*fx+fy*fy)
	xmn = min(x)
        xmx = max(x)
        ymn = min(y)
        ymax= max(y)
	xstep=(xmx-xmn)/50
	ystep=(ymax-ymn)/50
	
	h,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)])
	hvx,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=vx)
	hvy,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=vy)
	Jx =  hvx
	Jy =  hvy
	hvx = where(h>0,hvx/h,0)
	hvy = where(h>0,hvy/h,0)
	for i in arange (1,nf):					
		x,y=transpose(loadtxt("./posizione/dati_"+str(i)))
		vx,vy=transpose(loadtxt("./velocity/velocity_"+str(i)))
		h,bx_temp,by_temp = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)])
		hvx_temp,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=vx)
		hvy_temp,bx,by = histogram2d(x,y,bins=[arange(xmn,xmx,xstep),arange(ymn,ymax,ystep)],weights=vy)
		Jx  +=hvx_temp
		Jy  +=hvy_temp
		hvx_temp = where(h>0,hvx_temp/h,0)
		hvy_temp = where(h>0,hvy_temp/h,0)
		hvx+=hvx_temp
		hvy+=hvy_temp
	v=hvx**2+hvy**2
	
	#print(hvy)
	title("vector field")
	quiver(bx[:-1]+xstep,by[:-1]+ystep,hvx/float(nf),hvy/float(nf),v,scale=5)
	colorbar()
	figure()
	title("current field")
	quiver(bx[:-1]+xstep,by[:-1]+ystep,Jx/float(len(x)*nf),Jy/float(len(x)*nf),v,scale=0.002)
#        quiver(by[:-1],bx[:-1],hfy,hfx,f,scale=20)

	#print("la media del modulo della forza e' "+str(mean(modulo)))
	ax=gca()#fissa gli assi
	ax.set_aspect(1)#mette la stessa scala  tra x e y
	show()
	colorbar()	
	print("somma delle vy"+str(sum(vy)))
	print("somma delle vx"+str(sum(vx)))



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
	
		
