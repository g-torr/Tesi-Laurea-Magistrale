def stampa():
	b,sh=np.load("F_cumulative.npy")
	b,h=np.load("F.npy")
	D=0.1
	tau=0.06
	L=sqrt(D*tau)
	cond=(b>sin(pi/6)-1+L)&(b<-L)
	cond1=b>L
	p=polyfit(b[cond], sh[cond] ,1)
	M1=abs(min(h))
	m1=max((h))
	M=max(sh[b<0]-polyval(p,b)[b<0])
	p1=polyfit(b[cond1], sh[cond1] ,1)
	m=max(sh[b>0]-polyval(p1,b)[b>0])
	apex=2*D*L*cos(pi/6)
	base=2*D*L*2*cos(pi/4)*cos(pi/12)
	'''print(" forza alla base tramite" )
	print("cumulative		Fy			theoretical  ")
	print(str(m)+"		"+str(mean(h[cond]))+"		"+str(2*D*L*2*cos(pi/4)*cos(11*pi/12)))'''
	print(" forza al centro" )
	print("cumulative		Fy			teorica ")
	print(str(M)+"		"+str(max(h)-max(h[cond]))+"		"+str(2*D*L*cos(pi/3)))
	'''print("forza sulla punta")
	print("cumulative		Fy			teorica ")
	print(str(M)+"		"+str(min(h))+"		"+str(2*D*L*cos(pi/6)))'''
	plot(b,h)
	axhline(y=-2*D*0.001*(cos(pi/3)/sin(pi/3)-cos(pi/6)/sin(pi/6)),xmin=0.35,xmax=0.5,color="y")
	axhline(y=-2*D*0.001*(cos(pi/3)/sin(pi/3)),xmin=0.5,xmax=0.7,color="r")
	legend(["simulation","theory","theory"],"best")
	title ("$F_y$ outside")
	xlabel("$y$")
	ylabel("$F_y$")
