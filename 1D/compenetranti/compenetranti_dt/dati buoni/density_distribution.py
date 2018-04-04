
for dirname, dirnames, filenames in os.walk('.'):
	break
a=re.compile("dati")

for name in filenames:
	if a.match(name):	
		print("loading file"+ name)
		x=loadtxt(name)
		dt=float(re.sub("dati_","",name))
		h,b=histogram(x,bins=linspace(x.min(),x.max(),100))
		plot((b[1:]+b[:-1])/2.,h,"-.")
		legend(dt,loc="best")
