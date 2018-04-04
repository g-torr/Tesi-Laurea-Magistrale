from scipy.optimize import curve_fit
import scipy.special
def f(x,A,B):
	return A*x + B
def binna(r1,r2,D=0.2,tau=0.06,nb=10):
	x,y=np.load("probability_in_curv.npy")
	bns = linspace(x.min(),x.max(),nb)
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
	plot(hx,f(hx,A,B),"-",label="Best fit")
	m=D*tau/normalization(float(r1),float(r2),D,tau)
	q=sqrt(D*tau)/normalization(float(r1),float(r2),D,tau)
	plot(hx,hx*m+q,"-.",label="UCNA")
	legend(loc="best",numpoints=1,fontsize=12)#.draw_frame(False)
	xlabel("$\kappa$")
	ylabel("$P$")
	#title ("Probability vs curvature")
def normalization(r1,r2,D,tau):
	
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	return 	pi*r1*r2+sqrt(D*tau)*arc_lenght+D*tau*2*pi	
	
		
