def pldist(nome):
    dati = loadtxt(nome)
    h,b = histogram(dati, bins=linspace(dati.min(),dati.max(),50),normed=True)
    x = (b[1:]+b[:-1])/2. #credo serva per predendere il centro
    plot(x,h,'-ob')
    #-------------------------------------
    k = 2.
    tau = 0.6
    u = 0.5*k*x**2
    u1 = k*x
    D = 0.1
    y = exp(-u/D - tau*(u1**2)/(2.*D))
    y = y/trapz(y,x)
    plot(x,y,'-r')
    print "<x^2>="+str(trapz(y*x**2,x))

