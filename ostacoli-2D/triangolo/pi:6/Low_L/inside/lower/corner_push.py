D=0.1
L=sqrt(0.1*0.06)
b,sh=np.load("F_cumulative.npy")
theta=pi/6
cond1=((b-1/sin(theta)+1)<0)&((b-1/sin(theta)+1)>-1e-2)
mean(sh[cond1][:-2])-sh[-1]

