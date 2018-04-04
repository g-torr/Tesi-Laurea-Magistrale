def normalization(R=1.,theta=pi/6):
	D=0.1
	tau=0.06
	size=3.
	L=sqrt(D*tau)
	S=R**2*cos(theta)**3/sin(theta)-R**2/2*(pi-2*theta-sin(2*theta))
	Fy=2*D*(R+L)*cos(theta)-2*R*D*cos(theta)
	N=size**2-S+2*L*R/tan(theta)+L*R*(pi-2*theta)+L**2*(pi-2*theta)
	print("la UCNA ti da forza pari a "+str(Fy/N))
	print("frazione di particelle al bordo"+str((N-size**2+S)/N))
	#normal_force=
	#print("la forza normale ="+str(
	return  N	
	
