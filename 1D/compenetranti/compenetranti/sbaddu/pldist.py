def pldist(nome):
    x,v =transpose( loadtxt(nome))
    h,b = histogram(v, bins=linspace(v.min(),v.max(),50),normed=True)
    c= (b[1:]+b[:-1])/2. #credo serva per predendere il centro
    figure(1)
    title("velocity")
    plot(c,h,'-ob')
    h,b = histogram(x, bins=linspace(x.min(),x.max(),50),normed=True)
    c= (b[1:]+b[:-1])/2. #credo serva per predendere il centro
    figure(2)
    title("position")
    plot(c,h,'-ob')

    #-------------------------------------
'''    k = 2.
    tau = 0.6
    u = 0.5*k*x**2
    u1 = k*x
    D = 0.1
    y = exp(-u/D - tau*(u1**2)/(2.*D))
    y = y/trapz(y,x)
#    plot(x,y,'-r')
    print "<x^2>="+str(trapz(y*x**2,x))
'''
