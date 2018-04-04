def confronta():
	L=sqrt(0.5/6.)#la nuova lunghezza di persistenza e' v/D dell'active brownian in cui 
	angolo=[]
	value=[]
	D=0.5**2/(2*6.) #prendo la nuova diffusivita' pari a v**2/D basandomi solo sul calcolo dimensionale; considerando che la nuova diffusivita' e' sull'asse divido per 2
	figure(1,figsize=(4,3))

	'''
	b,sh=np.load("./pi:24/F_cumulative.npy")
	theta=pi/24
	angolo+=[theta]
	value=append(value,mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)	
	'''

	b,sh=np.load("./pi:12/F_cumulative.npy")
	theta=pi/12
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)	

	b,sh=np.load("./pi:6/F_cumulative.npy")
	theta=pi/6
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)	

	b,sh=np.load("./5pi:24/F_cumulative.npy")
	theta=5*pi/24
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)

	b,sh=np.load("./pi:4/F_cumulative.npy")
	theta=pi/4
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)
	
	b,sh=np.load("./7pi:24/F_cumulative.npy")
	theta=7*pi/24
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)

	b,sh=np.load("./pi:3/F_cumulative.npy")
	theta=pi/3
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)

	b,sh=np.load("./9pi:24/F_cumulative.npy")
	theta=9*pi/24
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-1]))
	plot(b,sh)

	b,sh=np.load("./5pi:12/F_cumulative.npy")
	theta=5*pi/12
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-2]))
	plot(b,sh)

	'''b,sh=np.load("./11pi:24/F_cumulative.npy")
	theta=11*pi/24
	cond=(b>L/2*cos(theta)+sin(theta))&(b<sin(theta)+L*cos(theta))
	cond1=b>sin(theta)+L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))'''

	legend(["$\\frac{\pi}{12}$","$\\frac{\pi}{6}$","$ \\frac{5\pi}{24}$","$\\frac{\pi}{4}$","$\\frac{7pi}{24}$","$\\frac{\pi}{3}$","$\\frac{3pi}{8}$","$\\frac{5pi}{12}$"],"best")
	title ("Cumulative sum  of the $F_y$")
	xlabel("h")
	ylabel("$\\sum_{j=0}^h F_y(j)$")
	#tight_layout()

	figure(2,figsize=(4,3))
	
	plot(cos(angolo),value,"-o")	
	t=linspace(0,1,100)
	plot(t,2*D*(L+1)*t)
	
	legend(["simulation","theoretical line"],"best")
	xlabel("$cos(\\theta)$")
	ylabel("Wedge push")
	title("Wedge push at different angles")
	tight_layout()

	np.save("confronta.npy",array([angolo,value]))

