def getlin(xl,yl,x,y):
	plot(xl,yl,"-")
	mx = xl[1]-xl[0]
	my = yl[1]-yl[0]
	tx = (x-xl[0])/mx
	ty = (y-yl[0])/my

