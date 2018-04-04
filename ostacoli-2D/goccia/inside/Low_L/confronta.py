def confronta():
	L=sqrt(0.1*0.06)
	angolo=[]
	value=[]


	b,sh=np.load("./theta=pi:12/point/F_cumulative.npy")
	theta=pi/12
	cond=(b>L+sin(theta))&(b<10*L+sin(theta))
	cond1=b>10*L+sin(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))

	b,sh=np.load("./theta=pi:8/point/F_cumulative.npy")
	theta=pi/8
	cond=(b>L+sin(theta))&(b<10*L+sin(theta))
	cond1=b>10*L+sin(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))

	b,sh=np.load("./theta=pi:5/point/F_cumulative.npy")
	theta=pi/5
	cond=(b>L+sin(theta))&(b<10*L+sin(theta))
	cond1=b>10*L+sin(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))

	b,sh=np.load("./theta=pi:4/point/F_cumulative.npy")
	theta=pi/4
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))

	b,sh=np.load("./theta=pi:3/point/F_cumulative.npy")
	theta=pi/3
	cond=(b>L+sin(theta))&(b<max(b)-2*L)
	cond1=b>max(b)-2*L
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))
	plot(cos(angolo),value,"-o")



	


	np.save("confronta.npy",array([angolo,value]))
