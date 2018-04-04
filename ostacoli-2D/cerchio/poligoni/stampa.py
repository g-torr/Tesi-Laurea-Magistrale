def	probability(D,tau):
	return 2*pi*sqrt(D*tau)*(1+sqrt(D*tau))/(pi+2*pi*sqrt(D*tau)*(1+sqrt(D*tau)))
def stampa():
	L,f1,p1=transpose(loadtxt("/home/giuseppe/Documents/myCUDA/ostacoli-2D/cerchio/poligoni/changing_L/sim_vs_L.dat"))
	tau,f,p=transpose(loadtxt("/home/giuseppe/Documents/myCUDA/ostacoli-2D/cerchio/poligoni/D=4/sim_vs_tau.dat"))
	t=linspace(0.05,16**2,1000)
	#m=probability(4,tau[1])/p[1]
	figure(figsize=(8,6))
	plot(t, probability(1,t),"b--",L**2,p1,"b.",t, probability(4,t),"m--",tau,p,"m.")
	legend(["UCNA Prediction $D=1$","Poligon Simulation $D=1$","UCNA Prediction $D=4$","Poligon Simulation $D=4$",],"best")
	semilogx()
	xlabel ("$\\tau$")
	ylabel ("$p$")
	title ("Probability density of being on the boundary vs $\\tau$")
