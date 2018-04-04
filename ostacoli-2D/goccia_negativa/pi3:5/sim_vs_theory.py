b,sh=np.load("F_cumulative.npy")
D=0.1
L=sqrt(0.1*0.06)
theta=pi/6
plot(b,sh,"b")
t=linspace(-1,0,100)
plot(t,-(sqrt(2*pi)*L+2*R)*D*cos(arcsin(t)),"g")
t=linspace(0,0.5,100)
F_negative_corner=sqrt(2*pi)*L*D*(cos(theta)-2*cos(pi/4))
plot(t,-(sqrt(2*pi)*L+2*R)*D*cos(arcsin(t))+ 2*D*R*t*cos(theta)/sin(theta)+F_negative_corner,"g")

