def confronta():
	L=sqrt(0.5/6.)
	angolo=[]
	value=[]


	b,h=np.load("./pi:12/F.npy")
	theta=pi/12
	'''cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L
	p=polyfit(b[cond], h[cond] ,1)'''
	angolo+=[theta]
	value=append(value,max(h))
	plot(b,h)	

	b,h=np.load("./pi:6/F.npy")
	theta=pi/6
	'''cond=(b>L+sin(theta))&(b<10*L+sin(theta))
	cond1=b>10*L+sin(theta)
	p=polyfit(b[cond], h[cond] ,1)'''
	angolo+=[theta]
	value=append(value,max(h))
	plot(b,h)

	b,h=np.load("./5pi:24/F.npy")
	theta=5*pi/24
	angolo+=[theta]
	value=append(value,max(h))
	plot(b,h)

	b,h=np.load("./pi:4/F.npy")
	theta=pi/4
	'''cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L
	p=polyfit(b[cond], h[cond] ,1)'''
	angolo+=[theta]
	value=append(value,max(h))
	plot(b,h)

	b,h=np.load("./7pi:24/F.npy")
	theta=7*pi/24
	angolo+=[theta]
	value=append(value,max(h))	
	plot(b,h)

	b,h=np.load("./pi:3/F.npy")
	theta=pi/3
	'''cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L
	p=polyfit(b[cond], h[cond] ,1)'''
	angolo+=[theta]
	value=append(value,max(h))	
	plot(b,h)

	b,h=np.load("./9pi:24/F.npy")
	theta=9*pi/24
	angolo+=[theta]
	value=append(value,max(h))	
	plot(b,h)
	'''b,h=np.load("./5pi:12/F.npy")
	theta=5*pi/12

	angolo+=[theta]
	value=append(value,max(h))	


	b,h=np.load("./11pi:24/F.npy")
	theta=11*pi/24
	angolo+=[theta]
	value=append(value,max(h))'''

	legend(["$\\frac{\pi}{12}$","$\\frac{\pi}{6}$","$ \\frac{5\pi}{24}$","$\\frac{\pi}{4}$","$\\frac{7pi}{24}$","$\\frac{\pi}{3}$","$\\frac{3pi}{8}$"],"best")
	figure()	


	plot(cos(angolo),value,"-o")
	D=0.1
	L=sqrt(D*0.06)
	t=linspace(0,1,100)
	plot(t,2*D*L*t)

	legend(["simulation","theoretical line"],"best")
	xlabel("$cos(\\theta)$")
	ylabel("Corner push")
	title("Corner push at different angles")


	


	np.save("confronta2.npy",array([angolo,value]))
