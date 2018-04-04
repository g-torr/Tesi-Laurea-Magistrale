def genera(theta=pi/6):
	x=[cos(theta),0,-cos(theta),0]
	y=[sin(theta)-1.,1/sin(theta)-1,sin(theta)-1.,0]
	plot(x,y)
	savetxt("input.dat",transpose(array([x,y])).astype(float32))
	print("il numero di vertici = "+str(len(x)))
