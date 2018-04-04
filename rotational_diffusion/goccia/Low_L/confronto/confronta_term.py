def confronta():
	D=0.5**2/(2.*6.) #prendo la nuova diffusivita' pari a v**2/D basandomi solo sul calcolo dimensionale; considerando che la nuova diffusivita' e' sull'asse divido per 2
	L=0.5/6.#la nuova lunghezza di persistenza e' v/D dell'active brownian in cui 
	angolo=[]
	value=[]
	figure(figsize=[4,3])

	'''
	b,sh=np.load("./pi:24/F_cumulative"+str(i)+".npy")
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
	sh_sum=0
	for i in arange(0,10):
		b,sh=np.load("./pi:12/F_cumulative"+str(i)+".npy")
		theta=pi/12
		F=-2*D*(1+L)*cos(theta)
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1/sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
		sh_sum=sh_sum+sh
	plot(b,sh/10.)

	angolo+=[theta]

	sh_sum=0
	for i in arange(0,10):
		b,sh=np.load("./pi:6/F_cumulative"+str(i)+".npy")
		theta=pi/6
		F=-2*D*(1+L)*cos(theta)
		correction=0
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1/sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
		sh_sum=sh_sum+sh
	plot(b,sh/10.)

	angolo+=[theta]

	sh_sum=0
	for i in arange(0,10):	
		b,sh=np.load("./5pi:24/F_cumulative"+str(i)+".npy")
		theta=5*pi/24
		F=-2*D*(1+L)*cos(theta)
		correction=0
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1/sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
		sh_sum=sh_sum+sh
	plot(b,sh/10.)

	angolo+=[theta]

	sh_sum=0
	for i in arange(0,10):	
		b,sh=np.load("./pi:4/F_cumulative"+str(i)+".npy")
		theta=pi/4
		F=-2*D*(1+L)*cos(theta)
		correction=0
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1/sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
		sh_sum=sh_sum+sh
	plot(b,sh/10.)

	angolo+=[theta]
	
	sh_sum=0
	for i in arange(0,10):	
		b,sh=np.load("./7pi:24/F_cumulative"+str(i)+".npy")
		theta=7*pi/24
		F=-2*D*(1+L)*cos(theta)
		correction=0
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1./sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
		sh_sum=sh_sum+sh
	plot(b,sh/10.)

	angolo+=[theta]

	sh_sum=0
	for i in arange(0,10):
		b,sh=np.load("./pi:3/F_cumulative"+str(i)+".npy")
		theta=pi/3
		F=-2*D*(1+L)*cos(theta)
		correction=0
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1/sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
		sh_sum=sh_sum+sh
	plot(b,sh/10.)

	angolo+=[theta]

	sh_sum=0
	for i in arange(0,10):
		b,sh=np.load("./9pi:24/F_cumulative"+str(i)+".npy")
		theta=9*pi/24
		F=-2*D*(1+L)*cos(theta)
		correction=0
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1/sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
		sh_sum=sh_sum+sh
	plot(b,sh/10.)
	angolo+=[theta]
	'''
	for i in arange(0,10):
		b,sh=np.load("./5pi:12/F_cumulative"+str(i)+".npy")
		theta=5*pi/12
		F=-2*D*(1+L)*cos(theta)
		correction=0
		cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
		cond1=b>1/sin(theta)-2*L*cos(theta)
		p=polyfit(b[cond], sh[cond] ,1)
		value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
	angolo+=[theta]

	b,sh=np.load("./11pi:24/F_cumulative"+str(i)+".npy")
	theta=11*pi/24
	cond=(b>L/2*cos(theta)+sin(theta))&(b<sin(theta)+L*cos(theta))
	cond1=b>sin(theta)+L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))'''

	#legend(["$\\frac{\pi}{12}$","$\\frac{\pi}{6}$","$ \\frac{5\pi}{24}$","$\\frac{\pi}{4}$","$\\frac{7pi}{24}$","$\\frac{\pi}{3}$","$\\frac{3pi}{8}$","$\\frac{5pi}{12}$"],"best")
	#plot(cos(angolo),value,"-o")
	t=linspace(0,1,100)
	#plot(t,2*D*L*t/1.5933035442369659)
	figure(figsize=[4,3])
	a=reshape(value,[7,10])
	media=mean(a,1)
	err=std(a,1)/sqrt(10)
	plot(cos(angolo),media/(2*D*L),"o",ms=6,mec="black",mew=2,c="w",label="Simulation")
	''' questi plot danno il mapping
	plot(t,t,"r--",label="Prediction")
	#plot(t,t*1.4523,"m:",label="Rescaling")
	plotto invece il fit con intercetta 0'''
	x = cos(angolo)[:,np.newaxis]
	a, _, _, _ = np.linalg.lstsq(x,media/(2*D*L) )
	plot(t,a*t,label="Linear trend")
	xlabel("$cos(\\theta)$",fontsize=8)
	ylabel("$F_{corner}/\\left(2D\\mathscr{L}\\rho_{0}\\right)$",fontsize=8)
	tight_layout()
	xlabel("$cos(\\theta)$",fontsize=16)
	ylabel("$F_{corner}/\\left(2D\\mathscr{L}\\rho_{B}\\right)$",fontsize=16)

	np.save("confronta.npy",array([angolo,value]))
	fill_between(cos(angolo),(media-err)/(2*D*L),(media+err)/(2*D*L),color="gray",alpha=0.5)
	legend(loc="best",numpoints=1, frameon=False)
	return cos(angolo),media/(2*D*L)

