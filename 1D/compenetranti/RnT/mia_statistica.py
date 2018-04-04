def istogrammi(nome, num_bins):
	dati=loadtxt(nome)
	h,b=histogram(dati,bins=linspace(dati.min()-0.1,dati.max()+0.1,num_bins),normed=True)
	x = (b[1:]+b[:-1])/2.
   	plot(x,h,'-ob')

def pallette_al_bordo(nome): # guardiamo alla frazione di pallette al bordo rispetto a quelle totali

	dati=loadtxt(nome)
	a=len(dati)
	print("il mio vettore ha"+str(a))
	destra=len(dati[dati==dati.max()])	
	print("a destra ne abbiamo"+str(destra))
	sinistra=len(dati[dati==dati.min()])
	print("a sinistra ne ho"+str(sinistra))
	print("il rapporto pallette ai bordi rispetto alle pallette totali="+str((float)(destra+sinistra)/a))
