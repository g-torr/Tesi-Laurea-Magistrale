D=0.1
tau=0.06
def stampa():
	lista=[]
	Raggi=[0.5,1,1.5,2,3,4,6,8,12]
	for R in Raggi:
		try:	
			R,p,std_p=np.load("./R="+str(R)+"/pressione.npy")
			lista+=[R,p,std_p]	
		except ValueError:
			R,p=np.load("./R="+str(R)+"/pressione.npy")
			print("R="+str(R)+"does not contain the variance")

	R,p,std_p=transpose(reshape(lista,[9,3]))
	plot(1/R,p/D-1,'o',ms=10,mec="k",mew=3,c="w")
	fit=polyfit(1/R,p/D-1., 1)
	plot(1/R,polyval(fit,1/R),label="fit")
	fill_between(1/R,(p-std_p)/D-1,(p+std_p)/D-1,color="gray",alpha=0.3)
	xlabel ("$\\frac{1}{R}$")
	ylabel("Pressure")
	r=linspace(0,1,100)
	plot(r,sqrt(pi/2)*sqrt(0.1*0.06)*r,"--",label="Theory")
	return R,p,std_p
	
