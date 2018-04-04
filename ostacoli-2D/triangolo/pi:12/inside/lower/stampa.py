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
	print("simulated force on apex througth" )
	print("cumulative		Fy			theoretical on apex ")
	print(str(m)+"		"+str(m1)+"		"+str(apex))
	print("cumulative		Fy			theoretical on apex ")
	print(str(M)+"		"+str(M1)+"		"+str(base))
	print("ratio	cumulative	Fy		theoretical")
	print("	"+str(M/m)+"	"+str(M1/m1)+"	"+str(base/apex))	

