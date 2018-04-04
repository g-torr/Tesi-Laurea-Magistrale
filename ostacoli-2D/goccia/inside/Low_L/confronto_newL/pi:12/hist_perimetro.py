def hist_curvilineo(lmin):
 x0,y0=transpose(loadtxt("input.dat"))
 x,y=transpose(loadtxt("./posizione/dati_9"))
 fx,fy=transpose(loadtxt("./forza/forza_9"))
 x0=append(x0,x0[0])
 y0=append(y0,y0[0])
 s=[]
 d=0
 x=x[(fx!=0)|(fy!=0)]
 y=y[(fx!=0)|(fy!=0)]
 '''x = 3.*(rand(int(100000))-0.5)
 y = 3.*(rand(int(100000))-0.5)'''
 for i in arange(1, len(x0)):

	'''t = ((P.x-vertice[0].x)*(vertice[1].x-vertice[0].x)+(P.y-vertice[0].y)*(vertice[0+1].y-vertice[0].y))/
								((vertice[0+1].x-vertice[0].x)*(vertice[0+1].x-vertice[0].x)+(vertice[0+1].y-vertice[0].y)*(vertice[0+1].y-vertice[0].y));

	'''


  	t=((x-x0[i-1])*(x0[i]-x0[i-1])+(y-y0[i-1])*(y0[i]-y0[i-1]))/((x0[i]-x0[i-1])**2+(y0[i]-y0[i-1])**2)   
 	#t = where(t>1,1.,where(t<0,0.,t))
	xl=(x0[i]-x0[i-1])*t+x0[i-1]
	yl=(y0[i]-y0[i-1])*t+y0[i-1]
	l=sqrt((y-yl)**2+(x-xl)**2)	
	cond=(t>0)&(t<1)&(l<lmin)
	#cond = l<lmin
	print("prima della selezione ho "+str(len(x))+"punti")
	xc=x[cond]
	yc=y[cond]
	print("dopo la selezione ho "+str(len(xc))+"punti")
	figure(1)
	plot(xc,yc,"o")
	plot([x0[i-1],x0[i]],[y0[i-1],y0[i]],"-.")
	s_temp=sqrt((yc-y0[i-1])**2+(xc-x0[i-1])**2)
	s=append(s,d + s_temp)
	d+=sqrt((y0[i]-y0[i-1])**2+(x0[i]-x0[i-1])**2)
	figure(0)
 h,b=histogram(s,200)
 b=(b[1:]+b[:-1])/2
 plot(b[1:-1],h[1:-1],"o")
 return s

'''	x0,y0=transpose(loadtxt("input.dat"))
	x,y=transpose(loadtxt("./posizione/dati_9"))
	fx,fy=transpose(loadtxt("./forza/forza_9"))
	x0=append(x0,x0[0])
	y0=append(y0,y0[0])
	l=[]
	x=x[(fx!=0)|(fy!=0)]
	y=y[(fx!=0)|(fy!=0)]
	figure(3)
	plot(x,y,".")
	d=0
	for i in arange(1, len(x0)):
		m=(y0[i]-y0[i-1])/(x0[i]-x0[i-1])
		cond= abs(y-(y0[i-1]+m*(x-x0[i-1])))< 0.01	
		yc=y[cond]
		xc=x[cond]
		d_temp=sqrt((y0[i]-y0[i-1])**2+(x0[i]-x0[i-1])**2)
		l_temp=sqrt((yc-y0[i-1])**2+(xc-x0[i-1])**2)
		figure(1)
		plot(xc[l_temp<d_temp],yc[l_temp<d_temp],".")		
		l_temp=l_temp[l_temp<d_temp]
		
		l=append(l,d+l_temp)
		d+=d_temp
		
		
	figure(0)
	h,b=histogram(l,100)
	b=(b[1:]+b[:-1])/2
	plot(b[1:-1],h[1:-1],"o")
	return l'''
