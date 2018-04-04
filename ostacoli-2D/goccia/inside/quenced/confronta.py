def confronta():
	L=sqrt(0.1*0.06)
	angolo=[]
	value=[]
	D=0.1
	tau=0.06
	v=sqrt(D/tau)
	'''
	b,sh=np.load("./pi:24/F_cumulative.npy")
	theta=pi/24
	F=-2*D*(1+L)*cos(theta)
	correction=0
	cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
	cond1=b>1/sin(theta)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
	plot(b,sh)	
	'''

	b,sh=np.load("./pi:12/F_cumulative.npy")
	theta=pi/12
	cond=b>sin(theta)
	cond1=b<1/sin(theta)
	M=max(sh[cond])
	m=mean(sh[cond&cond1])
	angolo+=[theta]
	value=append(value,M-m)
	plot(b,sh)	

	b,sh=np.load("./pi:8/F_cumulative.npy")
	theta=pi/8
	cond=b>sin(theta)
	cond1=b<1/sin(theta)
	M=max(sh[cond])
	m=mean(sh[cond&cond1])
	angolo+=[theta]
	value=append(value,M-m)
	plot(b,sh)	


	b,sh=np.load("./pi:6/F_cumulative.npy")
	theta=pi/6
	cond=b>sin(theta)
	cond1=b<1/sin(theta)
	M=max(sh[cond])
	m=mean(sh[cond&cond1])
	angolo+=[theta]
	value=append(value,M-m)
	plot(b,sh)	

	b,sh=np.load("./pi:4/F_cumulative.npy")
	theta=pi/4
	cond=b>sin(theta)
	cond1=b<1/sin(theta)
	M=max(sh[cond])
	m=mean(sh[cond&cond1])
	angolo+=[theta]
	value=append(value,M-m)
	plot(b,sh)	

	b,sh=np.load("./pi:3/F_cumulative.npy")
	theta=pi/3
	cond=b>sin(theta)
	cond1=b<1/sin(theta)
	M=max(sh[cond])
	m=mean(sh[cond&cond1])
	angolo+=[theta]
	value=append(value,M-m)
	plot(b,sh)	


	legend(["$\\frac{\pi}{12}$","$\\frac{\pi}{8}$","$ \\frac{\pi}{6}$","$\\frac{\pi}{4}$","$\\frac{pi}{3}$"],"best")
	figure(figsize=[4,3])	
	plot(cos(angolo),value,"o",ms=10,mec="black",mew=3,c="w",label="Simulation")
	t=linspace(0,1,100)
	plot(t,t*v/sqrt(2*pi),"r--",label="Theory")
	legend(loc="best",numpoints=1,frameon=False)
	xlabel("$cos(\\theta)$")
	ylabel("$F_{Corner}$")
	tight_layout()
	xlabel("$cos(\\theta)$",fontsize=16)
	ylabel("$F_{Corner}$",fontsize=16)

	#title("Corner push at different angles")
	np.save("confronta.npy",array([angolo,value]))

