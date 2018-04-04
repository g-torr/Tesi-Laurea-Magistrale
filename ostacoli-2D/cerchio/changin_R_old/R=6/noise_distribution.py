def noise_distribution(start,stop):
	x,y=transpose(loadtxt("./posizione/dati_"+str(start)))
	vx,vy=transpose(loadtxt("./velocity/velocity_"+str(start)))	
	fx,fy=transpose(loadtxt("./forza/forza_"+str(start)))	
	for i in arange(start+1,stop+1):
		x_temp,y_temp=transpose(loadtxt("./posizione/dati_"+str(i)))
		vx_temp,vy_temp=transpose(loadtxt("./velocity/velocity_"+str(i)))	
		fx_temp,fy_temp=transpose(loadtxt("./forza/forza_"+str(i)))	
		x=append(x,x_temp)
		y=append(y,y_temp)
		vx=append(vx,vx_temp)
		vy=append(vy,vy_temp)
		fx=append(fx,fx_temp)
		fy=append(fy,fy_temp)
	F=sqrt(fx**2+fy**2)
	cond= F>0
	eta_x=vx+fx
	eta_y=vy+fy
	eta_phi=(eta_x*y-eta_y*x)/(sqrt(x**2+y**2))
	eta_r=(eta_x*x+eta_y*y)/(sqrt(x**2+y**2))
	h,bx,by=histogram2d(eta_r[cond], eta_phi[cond],bins=50,normed="True")
	imshow(h,origin="low",extent=[by[0],by[-1],bx[0],bx[-1]])
	colorbar()
	figure()
	plot(eta_phi[cond], eta_r[cond], ".", alpha=0.01)
	xlabel("$\eta_{\phi}$")
	ylabel("$\eta_{r}$")
	figure()
	costheta=(eta_x*x+eta_y*y)/(sqrt(x**2+y**2)*sqrt(eta_x**2+eta_y**2))
	theta=arccos(costheta)
	h,b=histogram(theta[cond],32,normed="True")
	plot((b[:-1]+b[1:])/2,h)	
	xlabel("$\\theta$")
	ylabel("istogramma angolo compreso")
	

