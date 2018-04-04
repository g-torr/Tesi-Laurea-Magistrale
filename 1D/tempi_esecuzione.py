def tempi_esecuzione():
	figure(figsize=[4,3])
	dt=[0.4,0.2,0.1,0.05,0.025]
	t_h=[20437,36933,69502,134876,265846]
	t_h=np.array(t_h)
	plot(dt,t_h/1000.,"o",ms=6,mec="b",mew=2,c="w",label="Hard Pot.") 
	t_s=[27284,49118,93535,181681,358777]
	t_s=np.array(t_s)
	plot(dt,t_s/1000.,"s",ms=6,mec="g",mew=2,c="w",label="Soft Pot.") 
	xlabel("$\\Delta t /\\tau$",fontsize=16)
	ylabel("T.o.e",fontsize=16)
	tight_layout()
