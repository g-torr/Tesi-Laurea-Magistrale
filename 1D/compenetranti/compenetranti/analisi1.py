def analisi(fonam,nf):
  ms = []
  for i in arange(0,nf+1):
    fnam = fonam+"/forza_"+str(i)+".txt"
    data = loadtxt(fnam)
    datap = data[data>0]
    m = mean(datap)
    #print datap
    ms += [m]
  plot(ms)
  ms = array(ms).mean()
  return ms


def analisi_all(inizio=.5,n=11,nf=398): #nf number of files n = il numero di dt che gli passo
  dts =[]
  dtsf = []
  mss = []
  for i in arange(0,n):
	dt = inizio/2**i
	dtsf += [dt]
   	dt= '{0:.6f}'.format(dt)
   	dts += [dt]	
	#mss = []
	#for dt in dts:
	fonam = "forza"+str(dt)
	ms = analisi(fonam,nf) 
	mss += [ms] #append gli elementi alla lista
  mss = array(mss)
  dts = array(dtsf)
  title("simulazione con float")
  xlabel('dt')
  ylabel('forza media alle pareti')
  #gca().invert_xaxis()
  #loglog(dts,mss,"-o")
  semilogx(dts,mss,"-o")
  #return dts
