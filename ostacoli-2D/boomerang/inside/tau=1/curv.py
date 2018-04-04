from matplotlib.path import Path

def poly(ts,xs):
    ns = arange(0,len(xs))
    ys = []
    for t in ts:
        ys += [dot(xs,t**ns)]
    return array(ys)

def poly1(ts,xs):
    ns = arange(0,len(xs))
    ys = []
    for t in ts:
        myt = where(ns<=0,0,t**(ns-1))
        ys += [dot(xs,ns*myt)]
    return array(ys)

def poly2(ts,xs):
    ns = arange(0,len(xs))
    ys = []
    for t in ts:
        myt = where(ns-1<=0,0,t**(ns-2))
        ys += [dot(xs,(ns-1)*ns*myt)]
    return array(ys)

def curv(x1,y1,x2,y2):
    k = (x1*y2-y1*x2)/((x1**2+y1**2)**1.5)
    return k

def plcurve(D,tau,nop):
    xc = loadtxt("xc.dat")*10
    yc = loadtxt("yc.dat")*10
    ts = linspace(0,1,nop)
    xs = poly(ts,xc)
    ys = poly(ts,yc)
    xs1 = poly1(ts,xc)
    ys1 = poly1(ts,yc)
    xs2 = poly2(ts,xc)
    ys2 = poly2(ts,yc)
    figure(1)
    subplot(2,1,1)
    plot(ts,xs,'-')
    plot(ts,xs1,'--')
    plot(ts,xs2,'-.')
    subplot(2,1,2)
    plot(ts,ys,'-')
    plot(ts,ys1,'--')
    plot(ts,ys2,'-.')
    k = curv(xs1,ys1,xs2,ys2)
    figure(2)
    plot(ts,k)
    kmin = k.min()
    kmax = k.max()
    print "kmax = "+str(kmax)
    print "kmin = "+str(kmin)
    figure(3)
    for i in arange(0,len(xs)-1):
        myk = (k[i]+k[i+1])/2.
        myi = (myk-kmin)/(kmax-kmin)
        myc = (myi,0,1.-myi)
        myc = cm.jet(myi)
        plot(xs[i:i+2],ys[i:i+2],'-',lw=7,color=myc)
    cmx = 1.1*concatenate((abs(xs),abs(ys))).max()
    xlim(-cmx,cmx)
    ylim(-cmx,cmx)

import scipy.interpolate as scint
def 	xy(a,b,phi):
  	return a*np.cos(phi), b*np.sin(phi)

def splcurve(a,b):
	phis=np.arange(0,2*pi+pi/4,pi/4)
	x,y=xy(a,b,phis)    #da qui inizia la parte che mi interessa
	tck,u = scint.splprep([x,y],s=0)
	t1=arange(0,2*pi,pi/1)
	out = scint.splev(t1,tck)
	return out
	#plot(x,y,"-o")
'''
    k = curv(x1,y1,x2,y2)    
    kmin = k.min()
    kmax = k.max()
    for i in arange(0,len(xs)-1):
        myk = (k[i]+k[i+1])/2.
        myi = (myk-kmin)/(kmax-kmin)
        myc = (myi,0,1.-myi)
        #myc = cm.jet(myi)
	myc = cm.seismic(myi)
        plot(xs[i:i+2],ys[i:i+2],'-',lw=7,color=myc)
    cmx = 1.1*concatenate((abs(xs),abs(ys))).max()
	plot(x1,y1,"-o")
   # xlim(-cmx,cmx)
    #ylim(-cmx,cmx)
	ax = gca()
	ax.set_aspect(1)


'''


'''
def abs(x):
	k=x-1
	y=1+k
	#y=where(k<0,1./(1.-k),y)
	y=where(k<0,exp(k),y)
	return y'''
from matplotlib.path import Path
def area(xx,yy):
	'''t=linspace(0,2*pi,1000)
	xx=cos(t)
	yy=sin(t)'''
	myp=Path(transpose([xx,yy]))
	xmax=max(xx)
	ymax=max(yy)
	xmin=min(xx)
	ymin=min(yy)
	size=max([xmax,ymax,abs(xmin),abs(ymin)])
	myx= linspace(-size,size,1000)
	X,Y=meshgrid(myx,myx)
	c=myp.contains_points(transpose(array([X.flatten(),Y.flatten()])))
	#imshow(c.reshape(shape(X)))
	dxdy=diff(myx)[0]**2
	return sum(c)*dxdy
def boomerang(theta=pi/4):#out[0] ed out[1] sono rispettivamente ascisse e ordinate della mia curva
 x=[cos(theta),0,-cos(theta),0,cos(theta)]
 y=[sin(theta)-1.,1/sin(theta)-1,sin(theta)-1.,0,sin(theta)-1.]
 tck, u = scint.splprep([x, y],s=0,per=1)
 tnew = np.arange(0, 1.01, 0.01)
 out = scint.splev(tnew,tck)
 d1 = scint.splev(tnew,tck,der=1)
 d2 = scint.splev(tnew,tck,der=2)
 '''figure()
 plot(x, y, 'x', out[0], out[1], a*np.sin(2*np.pi*unew), b*np.cos(2*np.pi*unew), x, y, 'b')
 legend(['Linear', 'Cubic Spline', 'True'])
 axis([min(out[0])-0.5,max(out[0])+0.5,min(out[1])-0.5,max(out[1])+0.5 ])
 title('Spline of parametrically-defined curve')
 '''
      
 plot(out[0],out[1],x,y)
 #plot(tnew,k,tnew,a*b/(((a**2-b**2)*sin(2*pi*tnew)**2+b**2)**1.5))
 legend(['spline', 'teorica'])
 figure(2)  
 sign=input("press 1 if inside, -1 if outside")
 title("curvatura vs t")
 k = sign * curv(d1[0],d1[1],d2[0],d2[1])
 plot(tnew,k)
 somma=[[],[]]
 normal=[[],[]]
 tangenzial=[[],[]]
 figure(1)
 modulo=sqrt(d1[0]**2+d1[1]**2)
 dt=diff(tnew)[0]
 arc_lenght=sum(modulo*dt)
 superficie=area(out[0],out[1])
 D=.1
 tau=1.
 for i in arange (0,len(out[0])):  
 	#s=sqrt(out[0][i]**2+out[1][i]**2)
	tangenzial[0]=d1[0][i]/modulo[i]
	tangenzial[1]=d1[1][i]/modulo[i]
	normal[0]=np.array(normal[0])
	normal[1]=np.array(normal[1])
	normal[0]=sign*d1[1]/modulo
	normal[1]=-sign*d1[0]/modulo
	somma[0]=modulo*dt*normal[0]*D*abs(sqrt(D*tau)*k+1)
	somma[1]=modulo*dt*normal[1]*D*abs(sqrt(D*tau)*k+1)
	
 #return somma

 quiver(out[0],out[1],normal[0],normal[1])
 print("arc lenght="+str(arc_lenght)+" area= "+str(superficie))
 
 np.save("curvatura",k)
 kmin = k.min()
 kmax = k.max()
 figure(3)
 
 
 normalizzazione= superficie + sqrt(D*tau)*(abs(arc_lenght+sqrt(D*tau)*sum(k*modulo)*dt))
 print("UCNA theoretical results: \n normalization= "+str(normalizzazione))
 print("F_x risultante"+str(sum(somma[0])/normalizzazione)+"  F_y risultante"+str(sum(somma[1])/normalizzazione))
 particle_on_boundary=(normalizzazione-superficie)/normalizzazione
 print("  frazione di particelle sul bordo = "+str(particle_on_boundary))
 avg_module=abs(sum(modulo*dt*D*abs(sqrt(D*tau)*k+1)))/normalizzazione
 print(" average of the force module(respect to all particle) = "+str(avg_module))
 cond_avg_module=abs(sum(modulo*dt*D*abs(sqrt(D*tau)*k+1)))/(normalizzazione-superficie)
 print(" conditional average of the force module(respect to boundary) = "+str(cond_avg_module))
 for i in arange(0,len(out[0])-1):
        myk = (k[i]+k[i+1])/2.
        myi = (myk-kmin)/(kmax-kmin)
        myc = (myi,0,1.-myi)
        myc = cm.jet(myi)
	#myc = cm.seismic(myi)
        plot(out[0][i:i+2],out[1][i:i+2],'-',lw=7,color=myc)
 xmx = 1.1*max(out[0])
 ymx = 1.1*max(out[1])
 #plot(out[0],out[1],"-o")
 savetxt("input.dat",transpose(array([out[0][:-1],out[1][:-1]])).astype(float32))
 xlim(-xmx,xmx)
 ylim(-ymx,ymx)
 ax = gca()
 ax.set_aspect(1)
 figure(4)
 title("x(t),y(t)")
 plot(tnew,out[0],"-o") 
 plot(tnew,out[1],"-o")
