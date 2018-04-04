N=100000
def GCN(D=1,tau=1):
	figure()
	x=-1+2*rand(N)
	v=normal(0,D/tau,size=N)
	dt=0.1
	v=normal(0,1,size=N)
	for t in arange(0,10000):
	        v+=-v*dt/tau+normal(0,1,size=N)*sqrt(2*dt*D)/tau
	        x+=v*dt
	        x[x>1]=1
	        x[x<-1]=-1
	        #plot(t,x,"ko")    
	h,b=histogram(x,bins=linspace(x.min(),x.max(),100))
	b=(b[:-1]+b[1:])/2.
	xlim(0,)
	plot(b,h)
	title("GCN   D = "+str(D)+"  $\\tau = $"+str(tau))
def GrT(D=1,tau=1):
	figure()
	x=-1+2*rand(N)
	for t in arange(0,1000):
		v=normal(0,D/tau,size=N)
		x+=v*random.poisson(tau,size=N)
		x[x>1]=1
		x[x<-1]=-1
	h,b=histogram(x,bins=linspace(x.min(),x.max(),100))
	b=(b[:-1]+b[1:])/2.
	plot(b,h)
	#xlim(0,)
	title("Grt, D = "+str(D)+"  $\\tau = $"+str(tau))
def RnT():
	figure()
	x=-1+2*rand(N)
	tau=0.1	
	v=1
	for t in arange(0,1000):
		p=rand(N)
		x+=v*tau*sign(p-0.5)
		x[x>1]=1
		x[x<-1]=-1
	h,b=histogram(x,bins=linspace(x.min(),x.max(),100))
	b=(b[:-1]+b[1:])/2.
	plot(b,h)
	#xlim(0,)
	title("RnT")
