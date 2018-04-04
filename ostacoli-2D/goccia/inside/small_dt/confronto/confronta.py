def confronta():
	L=sqrt(0.1*0.06)
	angolo=[]
	value=[]


	b,sh=np.load("./pi:12/F_cumulative.npy")
	theta=pi/12
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(b,sh)	

	b,sh=np.load("./pi:6/F_cumulative.npy")
	theta=pi/6
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(b,sh)

	b,sh=np.load("./5pi:24/F_cumulative.npy")
	theta=5*pi/24
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(b,sh)

	b,sh=np.load("./pi:4/F_cumulative.npy")
	theta=pi/4
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(b,sh)
	
	b,sh=np.load("./7pi:24/F_cumulative.npy")
	theta=7*pi/24
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(b,sh)

	b,sh=np.load("./pi:3/F_cumulative.npy")
	theta=pi/3
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(b,sh)

	b,sh=np.load("./9pi:24/F_cumulative.npy")
	theta=9*pi/24
	cond=(b>sin(theta))&(b<2*L*cos(theta)+sin(theta))
	cond1=b>2*L*cos(theta)+sin(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(b,sh)

	'''b,sh=np.load("./5pi:12/F_cumulative.npy")
	theta=5*pi/12
	cond=(b>L/2*cos(theta)+sin(theta))&(b<sin(theta)+L*cos(theta))
	cond1=b>sin(theta)+L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))


	b,sh=np.load("./11pi:24/F_cumulative.npy")
	theta=11*pi/24
	cond=(b>L/2*cos(theta)+sin(theta))&(b<sin(theta)+L*cos(theta))
	cond1=b>sin(theta)+L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))'''

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

	np.save("confronta.npy",array([angolo,value]))

