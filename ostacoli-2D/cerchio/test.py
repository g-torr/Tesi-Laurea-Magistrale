from matplotlib.path import Path
import scipy.interpolate as scint
def curv(x1,y1,x2,y2):
    k = (x1*y2-y1*x2)/((x1**2+y1**2)**1.5)
    return k

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
def cerchio():#out[0] ed out[1] sono rispettivamente ascisse e ordinate della mia curva
 t=linspace(0,2*pi,100)
 x=cos(t)
 y=sin(t)
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
 tau=.06
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
 avg_module=abs(sum(modulo*dt*D*abs(sqrt(D*tau)*k+1)))/normalizzazione
 print(" average of the force module = "+str(avg_module))
 particle_on_boundary=(normalizzazione-superficie)/normalizzazione
 print("  frazione di particelle sul bordo = "+str(particle_on_boundary))
 
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
