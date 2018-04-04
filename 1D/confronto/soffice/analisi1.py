from pylab import*
def analisi(dt,nf):
   ms = []
   dt_s= '{0:.6f}'.format(dt)
   fonam = "forza"+dt_s
   devs=[]
   for i in arange(0,nf+1):
            fnam = fonam+"/forza_"+str(i)
            data = loadtxt(fnam)
            datap = abs(data)
            m = mean(datap)
            dev=var(datap)
    #print datap
            ms += [m]
            devs+=[dev]
 # plot(ms)
   ms = array(ms).mean()
   devs= sqrt(sum(devs))/(nf*len(devs))
   np.save(fonam, array([dt,ms]))
   np.save("errore_"+dt_s,array([dt,devs]))
   return ms


def analisi_all(inizio=.5,n=11,nf=497): #nf number of files n = il numero di dt che gli passo
  dts =[]
  dtsf = []
  mss = []
  for i in arange(0,n):
	dt = inizio/(1.2**i)
	dtsf += [dt]
   	#dt= '{0:.6f}'.format(dt)
   	#dts += [dt]	
	#mss = []
	#for dt in dts:
	ms = analisi(dt,nf) 
	mss += [ms] #append gli elementi alla lista
  mss = array(mss)
  dts = array(dtsf)
  title("simulazione con potenziale r 12")
  xlabel('dt')
  ylabel('forza media alle pareti')
  #gca().invert_xaxis()
  #loglog(dts,mss,"-o")
  semilogx(dts,mss,"-o")
  #mss=(mss-mss[-1])/mss[-1]
  np.save("from"+str(inizio)+"to"+str(dt)+".npy",array([dts,mss]))
