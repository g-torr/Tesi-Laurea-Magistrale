def genera(R,theta,N):
	t=linspace(pi-theta,theta,N)
	x=cos(t)
	y=sin(t)-R
	x=append(x,0)     #li usavo per poter fare il confronto con la goccia
	y=append(y,R/sin(theta)-R)
	#y=append(y,0.5)
	plot(x,y)	
	savetxt("input.dat",transpose(array([x,y])).astype(float32))
	print("il numero di vertici = "+str(len(x)))
