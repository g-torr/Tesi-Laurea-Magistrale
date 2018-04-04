import scipy.special
from scipy.optimize import curve_fit
import itertools
marker	= itertools.cycle(('o', 's', '^', '3')) 
color	= itertools.cycle(((0.8,0,0), (0.8,0,0.8), (0,0.8,0.8), (0.1,0.5,0.)))
r1 = 3.;r2 = 1.
D=[0.8,0.2,0.6]
tau=[0.03,0.06,0.12]
figure(figsize=[4,3])
def f(x,A,B):
	return A*x + B



def force_in_curvature():
	dati_x=[];dati_y=[]
	for j in arange(0,len(D)):
		L=sqrt(D[j]*tau[j])
		x=[];y=[];fx=[];fy=[];
		for i in arange(0,10):
			x_temp,y_temp=transpose(loadtxt("/home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/changing_tau/tau="+str(tau[j])+"/D="+str(D[j])+"/posizione/dati_"+str(i)))
			fx_temp,fy_temp=transpose(loadtxt("/home/giuseppe/Documents/myCUDA/ostacoli-2D/ellise/inside/changing_tau/tau="+str(tau[j])+"/D="+str(D[j])+"/forza/forza_"+str(i)))
			x=append(x,x_temp)
			y=append(y,y_temp)
			fx=append(fx,fx_temp)
			fy=append(fy,fy_temp)
		F=sqrt(fx**2+fy**2)
		r=sqrt(x**2+y**2)
		rho=count_nonzero(r<0.5)/(pi*0.5**2)
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
		hn = h/ls #e' la somma delle forze in d-theta per unita' di lunghezza /tutte le particelle
		curvatura=r1*r2*(9*sin(b[:-1])**2+cos(b[:-1])**2)**(-1.5)	
		hx,hy=binna(L*curvatura,hn/D[j]/rho-1.,7,L)
		dati_x=append(dati_x,hx)
		dati_y=append(dati_y,hy)
	
	dati_x,dati_y=sort([dati_x,dati_y])	
	np.save("dati.npy",[dati_x,dati_y])	
	p=polyfit(dati_x,dati_y,1)
	t=linspace(0,3,1000)
	plot(t*L,polyval(p,t*L),"k-.",label="Fit")
	plot(t*L,L*t,"r--",label="SCA")
	xlabel("$\\mathscr{L}/R$",fontsize=9)
	ylabel("$(\\hat{p}-\\hat{p}_0)/\\hat{p}_0$",fontsize=9)
	legend(loc=(0.,0.32),numpoints=1,frameon=False)
	tight_layout()
	xlabel("$\\mathscr{L}/R$",fontsize=16)
	ylabel("$(\\hat{p}-\\hat{p}_0)\\hat{p}_0$",fontsize=9)
	xlim(0.01,1)
	ylim(0.025,2)
	return dati_x,dati_y

def binna(x,y,numb_bins,L):
	bns = linspace(x.min(),x.max(),numb_bins)
	h,b = histogram(x,bins=bns)
	hx,b = histogram(x,bins=bns,weights=x)
	hy,b = histogram(x,bins=bns,weights=y)
	hy2,b = histogram(x,bins=bns,weights=y**2)
	hy2=hy2/h
	hx=hx/h
	hy=hy/h
	hbar=sqrt(hy2-hy**2)/sqrt(h) #std	
	mec=color.next()
	#fill_between(hx,hy-hbar,hy+hbar,color=tuple(sum(((1,1,1),mec),0)/2.),alpha=0.6)
	plot(hx,hy,marker.next(),ms=5,mec=mec,mew=2,c="w",label="$\mathscr{L}=$"+str(round(L,2)))
	errorbar(hx,hy,hbar,ecolor=mec,fmt=None)
	'''p=polyfit(hx,hy,1)
	plot(hx,polyval(p,hx),"-",label="Fit")'''
	return hx,hy
	#print("fit con y = "+str(A)+"x+"+str(B))


def normalization(r1,r2,D,tau):
	arc_lenght=4*r1*scipy.special.ellipe((r1**2-r2**2)/r1**2)
	return 	pi*r1*r2+sqrt(D*tau)*arc_lenght+D*tau*2*pi	



