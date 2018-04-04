import scipy.special
D=0.1;tau=0.06
def plell():
	r1 = 3.
	r2 = 1.
	D=0.1;tau=0.06
	x=[];y=[];fx=[];fy=[];
	for i in arange(0,20):
		x_temp,y_temp=transpose(loadtxt("./posizione/dati_"+str(i)))
		fx_temp,fy_temp=transpose(loadtxt("./forza/forza_"+str(i)))
		x=append(x,x_temp)
		y=append(y,y_temp)
		fx=append(fx,fx_temp)
		fy=append(fy,fy_temp)
	cond=fx**2+fy**2>0
	bns=linspace(-pi,pi,97)	
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
	hn = h/ls/len(x)
	a=linspace(min(h),max(h),100)
	imshow(transpose([a]),cmap="jet",aspect="auto",origin='lower')
	plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
	colorbar()


	f,(ax1,ax2)=subplots(2,1,figsize=(4,3),sharex='all',sharey='all')
	for i in arange(0,len(h)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		myi = (hn[i]-hn.min())/(hn.max()-hn.min()).astype(float)
		myc = (myi,0,1-myi)	
		ax1.plot([x0,x1],[y0,y1],'-',lw=6,color=cm.jet(myi))
	xlim(-r1*1.1,r1*1.1)
	ylim(-1.1*r2,1.1*r2)
	ax = gca()
	ax1.set_aspect(1)

	L=sqrt(D*tau)
	hU=L+L**2/(r1**2*sin(b)**2+r2**2*cos(b)**2)**1.5
	for i in arange(0,len(h)):
		x0 = r1*cos(b[i])
		x1 = r1*cos(b[i+1])
		y0 = r2*sin(b[i])
		y1 = r2*sin(b[i+1])
		myi = ((hU[i]+hU[i+1])/2.-hU.min())/(hU.max()-hU.min()).astype(float)
		ax2.plot([x0,x1],[y0,y1],'-',lw=6,color=cm.jet(myi))
	xlim(-r1*1.1,r1*1.1)
	ylim(-1.1*r2,1.1*r2)
	ax = gca()
	ax2.set_aspect(1)
	xlabel("$x$")
	ylabel("$y$")
	tight_layout()


	figure(3)
	curvatura=r1*r2*(9*sin(b[:-1])**2+cos(b[:-1])**2)**(-1.5)
	plot(curvatura,hn,".",label="Simulation")
	t=linspace(0,1,100)
	plot(t,sqrt(D*tau)*(1+sqrt(D*tau)*t)/normalization(3.,1.,D,tau))

	return hn
def force_in_curvature():
	r1 = 3.
	r2 = 1.
	D=0.1;tau=0.06
	x=[];y=[];fx=[];fy=[];
	for i in arange(0,20):
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
	hn = h/ls/len(x)
	curvatura=r1*r2*(9*sin(b[:-1])**2+cos(b[:-1])**2)**(-1.5)
	plot(curvatura,hn,".",label="Simulation")
	t=linspace(0,1,100)
	plot(t,D*(1+sqrt(D*tau)*t)/normalization(3.,1.,D,tau))
	return hn

def normalization(r1,r2,D,tau):
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	return 	pi*r1*r2+sqrt(D*tau)*arc_lenght+D*tau*2*pi	
	

