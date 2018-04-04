def forza(start,finish):
	fx,fy=transpose(loadtxt("./forza/forza_0"))
	modulo =sqrt(fx**2+fy**2)
	for i in arange(start+1,finish+1):
		fx,fy=transpose(loadtxt("./forza/forza_"+str(i)))
		temp_mod=sqrt(fx**2+fy**2)
		F=append(modulo,temp_mod)
	mean(F[F!=0])


fx,fy=transpose(loadtxt("./forza/forza_9"))
F=sqrt(fx**2+fy**2)
np.save("F_per_part.npy",mean(F[F!=0]))

