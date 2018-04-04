def confronta():
	L=sqrt(0.1*0.06)
	angolo=[]
	value=[]
	D=0.1
	figure(1,figsize=(4,3))

	'''
	b,sh=np.load("./pi:24/F_cumulative.npy")
	theta=pi/24
	angolo+=[theta]
	value=append(value,mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)	
	'''

	b,sh=np.load("./pi:12/F_cumulative.npy")
	theta=pi/12
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)	

	b,sh=np.load("./pi:6/F_cumulative.npy")
	theta=pi/6
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)	

	b,sh=np.load("./5pi:24/F_cumulative.npy")
	theta=5*pi/24
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)

	b,sh=np.load("./pi:4/F_cumulative.npy")
	theta=pi/4
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)
	
	b,sh=np.load("./7pi:24/F_cumulative.npy")
	theta=7*pi/24
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)

	b,sh=np.load("./pi:3/F_cumulative.npy")
	theta=pi/3
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)

	b,sh=np.load("./9pi:24/F_cumulative.npy")
	theta=9*pi/24
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)

	'''b,sh=np.load("./5pi:12/F_cumulative.npy")
	theta=5*pi/12
	angolo+=[theta]
	value=append(value,-mean(sh[abs(b-sin(theta))<1e-3]))
	plot(b,sh)

	b,sh=np.load("./11pi:24/F_cumulative.npy")
	theta=11*pi/24
	cond=(b>L/2*cos(theta)+sin(theta))&(b<sin(theta)+L*cos(theta))
	cond1=b>sin(theta)+L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))'''

	legend(["$\\pi/{12}$","${\pi}/{6}$","$ {5\pi}/{24}$","${\pi}/{4}$","${7\pi}/{24}$","${\pi}/{3}$","${3\pi}/{8}$"],"best",frameon=False)
	xlabel("$h$",fontsize=7)
	ylabel("$ F_y(h)/(v \\rho_{Bulk})$",fontsize=7)
	tight_layout()
	xlabel("$h$",fontsize=16)
	ylabel("$ F_y(h)/(v \\rho_{Bulk})$",fontsize=16)
	figure(2,figsize=(4,3))
	
	plot(cos(angolo),value/(2.*D*(L+1.)),"o",ms=5,mec="b",mew=2,c="w",alpha=0.7,label="Simulation")
	t=linspace(0,1,100)
	plot(t,t,"r--",label="Theory")
	xlabel("$\cos(\\theta)$",fontsize=16)
	ylabel("$F_{wedge}/\\left(2D\\left[\\mathscr{L}+R\\right]\\rho_{bulk}\\right)$",fontsize=16)
	legend(loc="best",numpoints=1,frameon=False)
	ylim(0,1.05)
	tight_layout()
	return array([angolo,value])
	#np.save("confronta.npy",array([angolo,value]))

