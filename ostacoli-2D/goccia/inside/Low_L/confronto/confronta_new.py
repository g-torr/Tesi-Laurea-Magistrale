def confronta():
	L=sqrt(0.1*0.06)
	angolo=[]
	value=[]
	D=0.1
	r=0


	b,sh=np.load("./pi:12/F_cumulative.npy")
	theta=pi/12
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][:-2])-sh[-1])
	plot(b,sh)	
	angolo+=[theta]
	#cond2=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)



	b,sh=np.load("./pi:6/F_cumulative.npy")
	theta=pi/6
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][:-2])-sh[-1])
	plot(b,sh)
	angolo+=[theta]

	'''(b-(sin(theta)-1)>0)&(b-(sin(theta)-1))'''

	b,sh=np.load("./5pi:24/F_cumulative.npy")
	theta=5*pi/24
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][:-2])-sh[-1])
	plot(b,sh)
	angolo+=[theta]


	b,sh=np.load("./pi:4/F_cumulative.npy")
	theta=pi/4
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][:-2])-sh[-1])
	plot(b,sh)
	angolo+=[theta]

	
	b,sh=np.load("./7pi:24/F_cumulative.npy")
	theta=7*pi/24
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][:-2])-sh[-1])
	plot(b,sh)
	angolo+=[theta]

	b,sh=np.load("./pi:3/F_cumulative.npy")
	theta=pi/3
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][:-2])-sh[-1])
	plot(b,sh)
	angolo+=[theta]

	b,sh=np.load("./9pi:24/F_cumulative.npy")
	theta=9*pi/24
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][:-2])-sh[-1])
	plot(b,sh)
	angolo+=[theta]

	
	'''b,sh=np.load("./5pi:12/F_cumulative.npy")
	theta=5*pi/12
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1][2:-2]))
	plot(b,sh)
	angolo+=[theta]

	b,sh=np.load("./11pi:24/F_cumulative.npy")
	theta=11*pi/24
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=((b-1/sin(theta)+r*sin(theta))<0)&((b-1/sin(theta)+r*sin(theta))>-1e-2)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,-mean(sh[cond1][2:-2]))'''

	legend(["$\\frac{\pi}{12}$","$\\frac{\pi}{6}$","$ \\frac{5\pi}{24}$","$\\frac{\pi}{4}$","$\\frac{7pi}{24}$","$\\frac{\pi}{3}$","$\\frac{3pi}{8}$","$\\frac{5pi}{12}$"],"best")
	figure(figsize=(4,3))	
	plot(cos(angolo),value/(D*L),"o",ms=6,mec="black",mew=2,c="w",label="Simulation")
	D=0.1
	L=sqrt(D*0.06)
	t=linspace(0,1,100)

	plot(t,sqrt(2*pi)*t,label="Prediction")
	legend(loc="best",numpoints=1, frameon=False)
	xlabel("$cos(\\theta)$",fontsize=8)
	ylabel("$F_{corner}/\\left(D\\mathscr{L}\\rho_{bulk}\\right)$",fontsize=10)
	tight_layout()
	xlabel("$cos(\\theta)$",fontsize=15)
	ylabel("$F_y^{corner}/\\left(D\\mathscr{L}\\rho_{B}\\right)$",fontsize=15)

	np.save("confronta.npy",array([angolo,value]))

