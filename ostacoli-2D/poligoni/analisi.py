def	grafica(i=8,b="input"):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	x1,y1=transpose(loadtxt(b))
	x1=append(x1,x1[0])
	y1=append(y1,y1[0])
	plot(x,y,"bo")
	plot(x1,y1,"-go")
	#xlim(-2,2)	
	
def	grafica2(i=8,b="input"):
	posizione= "./posizione/dati_"+str(i)
	x,y=transpose(loadtxt(posizione))
	x1,y1=transpose(loadtxt(b))
	x1=append(x1,x1[0])
	y1=append(y1,y1[0])
	#plot(x,y,"bo")
	xmn = x.min()
	xmx = x.max()
	h,bx,by = histogram2d(x,y,bins=arange(xmn,xmx,0.1))	
	h=transpose(h)	
	imshow(h,extent=(xmn,xmx,xmn,xmx))
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
	D=0.1	
	tau=0.6
	fmax=4*sqrt(D/tau)
	h_sum = forza(0,fmax,nb)
	for i in arange(1,8):
		h=forza(i,fmax,nb)
		h_sum+=h
	#h_sum=array(h_sum)
	b=linspace(-fmax,fmax,nb)
	x = (b[1:]+b[:-1])/2. #credo serva per predendere il centro
	plot(x,h_sum,'-k',lw=2)
	print("len(h)="+str(len(h)))
	print("len(b)="+str(len(b)))
