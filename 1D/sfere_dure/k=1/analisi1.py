def analisi(dt,nf):
  ms = []
  for i in arange(0,nf+1):
    dt_s= '{0:.6f}'.format(dt)
    fonam = "forza"+dt_s
    fnam = fonam+"/forza_"+str(i)
    data = loadtxt(fnam)
    datap = abs(data)
    m = mean(datap)
    #print datap
    ms += [m]
 # plot(ms)
  ms = array(ms).mean()
  np.save(fonam, array([dt,ms]))
  return ms


def analisi_all(inizio=.5,n=11,nf=96): #nf number of files n = il numero di dt che gli passo
  dts =[]
  dtsf = []
  mss = []
  for i in arange(0,n):
	dt = inizio/2**i
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
  mss=(mss-mss[-1])/mss[-1]
  np.save("sfere.npy",array([dts,mss]))
