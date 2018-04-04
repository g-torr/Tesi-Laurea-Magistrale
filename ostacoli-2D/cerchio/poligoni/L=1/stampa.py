def stampa():
	figure(figsize=(8,6))
	t=linspace(0.2,64,100)
	plot(t,2*pi*(1+1)/(pi+2*pi*(1+1))*t/t,"g-",tau,p,"g.")
	semilogx()
	xlabel("$\\tau$")
	ylabel("$p$")
	title ("Probability of being on boundary vs $\\tau$ at fixed $\mathscr{L}$")
	legend(["UCNA Prediction"," Simulation"],"best")

