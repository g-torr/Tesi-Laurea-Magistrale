import os
import sys
import re

for dirname, dirnames, filenames in os.walk('.'):
    break
def confronta():
	L=sqrt(0.1*0.06)
	raggi=[]
	value=[]
	D=0.1
	theta=pi/4.
	for name in dirnames:
		if "r=" in name:
			r=float(name[2:])
			b,sh=np.load(name+"/F_cumulative.npy")	
			cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
			cond1=abs(b-((1.-r)/sin(theta)+r*sin(theta)))<1e-2
			#cond1=(b-((1.-r)/sin(theta)+r*sin(theta))<0)&(b-((1.-r)/sin(theta)+r*sin(theta))>-1e-1)

			p=polyfit(b[cond], sh[cond] ,1)
			value=append(value,-mean(sh[cond1]))
			raggi+=[r]
	b,sh=np.load("point/F_cumulative.npy")
	r=0
	b,sh=np.load(name+"/F_cumulative.npy")	
	cond=(b>L*cos(theta)+sin(theta))&(b<(1.-r)/sin(theta)+r*sin(theta))
	cond1=abs(b-((1.-r)/sin(theta)+r*sin(theta))-L*cos(theta))<1e-1
	p=polyfit(b[cond], sh[cond] ,1)
	value=append(value,-mean(sh[cond1]))
	raggi+=[r]


	'''b,sh=np.load("./5pi:12/F_cumulative.npy")
	theta=5*pi/12
	F=-2*D*(1+L)*cos(theta)
	correction=0
	cond=(b>L*cos(theta)+sin(theta))&(b<1/sin(theta)-2*L*cos(theta))
	cond1=b>1/sin(theta)-2*L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1])-correction)
	plot(b,sh)

	b,sh=np.load("./11pi:24/F_cumulative.npy")
	theta=11*pi/24
	cond=(b>L/2*cos(theta)+sin(theta))&(b<sin(theta)+L*cos(theta))
	cond1=b>sin(theta)+L*cos(theta)
	p=polyfit(b[cond], sh[cond] ,1)
	angolo+=[theta]
	value=append(value,max(sh[cond1]-polyval(p,b)[cond1]))

	legend(["$\\frac{\pi}{12}$","$\\frac{\pi}{6}$","$ \\frac{5\pi}{24}$","$\\frac{\pi}{4}$","$\\frac{7pi}{24}$","$\\frac{\pi}{3}$","$\\frac{3pi}{8}$"],"best")'''
	figure(figsize=[4,3])	
	plot(raggi,value,"o",ms=7,mec="k",mew=2,c="w",label="Simulation")
	D=0.1
	L=sqrt(D*0.06)
	t=linspace(0,0.5,100)
	#plot(t,2*D*L*t/1.5933035442369659)
	plot(t,(2*t+sqrt(2*pi)*L)*D*cos(pi/4),label="Theory")
	#p1=polyfit(raggi,value,1)
	#plot(raggi,polyval(p1,raggi),"--",label="Fit")
	legend(frameon=False,loc="best")
	xlabel("$r$")
	ylabel("cup push")
	tight_layout()
	xlabel("$r$",fontsize=16)
	ylabel("cup push",fontsize=16)
	
