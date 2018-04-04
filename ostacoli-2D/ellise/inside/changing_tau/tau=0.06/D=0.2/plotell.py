import scipy.special
from scipy.optimize import curve_fit
D=0.2;tau=0.06;	r1 = 3.;r2 = 1.
def f(x,A,B):
	return A*x + B

def plell():
	x=[];y=[];fx=[];fy=[];
	for i in arange(0,10):
		x_temp,y_temp=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx_temp,fy_temp=transpose(loadtxt("./forza/forza_"+str(i)))
		x=append(x,x_temp)
		y=append(y,y_temp)
		fx=append(fx,fx_temp)
		fy=append(fy,fy_temp)
	cond=fx**2+fy**2>0
	bns=linspace(-pi,pi,96)	
	theta=arctan2(y[cond]/r2,x[cond]/r1)
	h,b=histogram(theta,bins=bns)
	ls = []
	for i in arange(0,len(h)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		#l = r1*abs(scipy.special.ellipeinc(arctan(tan(b[i+1]/3)),(r1**2-r2**2)/r1**2)-scipy.special.ellipeinc(arctan(tan(b[i]/3)),(r1**2-r2**2)/r1**2))# se volessi fare le cose per bene dovrei considerare la lunghezza ellittica, il problema e' che l'angolo da mettere non e' quello della parametrizzazione ma arctan(tan)b[i])), ma ho il problema del quadrante
		l = sqrt((x1-x0)**2+(y1-y0)**2)
		ls += [l]
	ls = array(ls)
	hn = h/ls/float(len(x[cond]))
	plotta(b,hn)
	figure()
	curvatura=r1*r2*(9*sin(b[:-1])**2+cos(b[:-1])**2)**(-1.5)
	binna(curvatura,hn)
	t=linspace(0,3,100)
	plot(t,(sqrt(D*tau)+(D*tau)*t)/(normalization(3.,1.,D,tau)-pi*r1*r2),"r-.",label="SCA")
	legend(loc=2,numpoints=1,frameon=False)
	xlabel("Curvature")
	ylabel("Surface Density")
	ylim(0.07,0.096)

	tight_layout()

	return hn
def force_in_curvature():
	x=[];y=[];fx=[];fy=[];
	for i in arange(0,10):
		x_temp,y_temp=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx_temp,fy_temp=transpose(loadtxt("./forza/forza_"+str(i)))
		x=append(x,x_temp)
		y=append(y,y_temp)
		fx=append(fx,fx_temp)
		fy=append(fy,fy_temp)
	F=sqrt(fx**2+fy**2)
	cond=F>0
	bns=linspace(-pi,pi,97)	
	theta=arctan2(y[cond]/r2,x[cond]/r1)
	h,b=histogram(theta,bins=bns,weights=F[cond])
	ls = []
	for i in arange(0,len(h)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		#l = r1*abs(scipy.special.ellipeinc(arctan(tan(b[i+1]/3)),(r1**2-r2**2)/r1**2)-scipy.special.ellipeinc(arctan(tan(b[i]/3)),(r1**2-r2**2)/r1**2))# se volessi fare le cose per bene dovrei considerare la lunghezza ellittica, il problema e' che l'angolo da mettere non e' quello della parametrizzazione ma arctan(tan)b[i])), ma ho il problema del quadrante
		l = sqrt((x1-x0)**2+(y1-y0)**2)
		ls += [l]
	ls = array(ls)
	hn = h/ls/float(len(x[cond]))
	plotta_forze(b,hn)
	curvatura=r1*r2*(9*sin(b[:-1])**2+cos(b[:-1])**2)**(-1.5)
	binna(curvatura,hn)
	t=linspace(0,3,100)
	plot(t,D*(1+sqrt(D*tau)*t)/(normalization(3.,1.,D,tau)-pi*r1*r2),"r-.",label="SCA")
	xlabel("Curvature")
	ylabel("Pressure")
	legend(loc=2,numpoints=1,frameon=False)
	ylim(0.127,0.192)
	tight_layout()

	return hn
def binna(x,y):
	bns = linspace(x.min(),x.max(),7)
	h,b = histogram(x,bins=bns)
	hx,b = histogram(x,bins=bns,weights=x)
	hy,b = histogram(x,bins=bns,weights=y)
	hy2,b = histogram(x,bins=bns,weights=y**2)
	hy2=hy2/h
	hx=hx/h
	hy=hy/h
	hbar=sqrt(hy2-hy**2)/sqrt(h) #std	
	plot(x,y,".")
	figure(figsize=[4,3])
	fill_between(hx,hy-hbar,hy+hbar,color=(.7,.7,.7))
	plot(hx,hy,"o",ms=5,mec="b",mew=2,c="w",alpha=0.7,label="Simulation")
	A,B= curve_fit(f, hx,hy)[0]
	plot(hx,f(hx,A,B),"k--",label="Fit")
	print("fit con y = "+str(A)+"x+"+str(B))

def plotta(b,hn):

	a=linspace(min(hn),max(hn),100)
	'''imshow(transpose([a]),cmap="jet",aspect="auto",origin='lower')
	plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
	colorbar()'''


	f,(ax1,ax2)=subplots(2,1,figsize=(4,3),sharex='all',sharey='all')
	M=hn.max();m=hn.min()
	for i in arange(0,len(hn)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		myi = (hn[i]-m)/(M-m).astype(float)
		myc = (myi,0,1-myi)	
		ax1.plot([x0,x1],[y0,y1],'-',lw=6,color=cm.jet(myi))
	xlim(-r1*1.1,r1*1.1)
	ylim(-1.1*r2,1.1*r2)
	ax = gca()
	ax1.set_aspect(1)
	L=sqrt(D*tau)
	hU=(L+L**2*r1*r2/(r1**2*sin(b)**2+r2**2*cos(b)**2)**1.5)/(normalization(3.,1.,D,tau)-pi*r1*r2)
	for i in arange(0,len(hn)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		myi = (hU[i]-m)/(M-m).astype(float)
		ax2.plot([x0,x1],[y0,y1],'-',lw=6,color=cm.jet(myi))
	xlim(-r1*1.1,r1*1.1)
	ylim(-1.1*r2,1.1*r2)
	ax = gca()
	ax2.set_aspect(1)
	xlabel("$x$")
	ylabel("$y$")
	tight_layout()
	xlabel("$x$",fontsize=16)
	ylabel("$y$",fontsize=16)
	f.subplots_adjust(right=0.75)
	ax1 = f.add_axes([0.77, 0.2, 0.05, 0.75])
	norm =Normalize(m,M)
	matplotlib.colorbar.ColorbarBase(ax1, cmap='jet',norm=norm,orientation='vertical')


def plotta_forze(b,hn):

	a=linspace(min(hn),max(hn),100)
	'''imshow(transpose([a]),cmap="jet",aspect="auto",origin='lower')
	plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
	colorbar()'''


	f,(ax1,ax2)=subplots(2,1,figsize=(4,3),sharex='all',sharey='all')
	M=hn.max();m=hn.min()
	for i in arange(0,len(hn)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		myi = (hn[i]-m)/(M-m).astype(float)
		myc = (myi,0,1-myi)	
		ax1.plot([x0,x1],[y0,y1],'-',lw=6,color=cm.jet(myi))
	xlim(-r1*1.1,r1*1.1)
	ylim(-1.1*r2,1.1*r2)
	ax = gca()
	ax1.set_aspect(1)

	L=sqrt(D*tau)
	hU=D*(1+L*r1*r2/(r1**2*sin(b)**2+r2**2*cos(b)**2)**1.5)/(normalization(3.,1.,D,tau)-pi*r1*r2)
	for i in arange(0,len(hn)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		myi = (hU[i]-m)/(M-m).astype(float)
		ax2.plot([x0,x1],[y0,y1],'-',lw=6,color=cm.jet(myi))
	xlim(-r1*1.1,r1*1.1)
	ylim(-1.1*r2,1.1*r2)
	ax = gca()
	ax2.set_aspect(1)
	xlabel("$x$")
	ylabel("$y$")
	tight_layout()
	xlabel("$x$",fontsize=16)
	ylabel("$y$",fontsize=16)

	f.subplots_adjust(right=0.75)
	ax1 = f.add_axes([0.77, 0.2, 0.05, 0.75])
	norm =Normalize(m,M)
	matplotlib.colorbar.ColorbarBase(ax1, cmap='jet',norm=norm,orientation='vertical')

def normalization(r1,r2,D,tau):
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	return 	pi*r1*r2+sqrt(D*tau)*arc_lenght+D*tau*2*pi	



