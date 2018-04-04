from pylab import*
from genera_goccia import*
def shapes():
	genera(1,0.7,pi/4)
	genera(1,0.6,pi/4)
	genera(1,0.5,pi/4)
	genera(1,0.4,pi/4)
	genera(1,0.3,pi/4)
	legend(["$r=0.7$","$r=0.6$","$r=0.5$","$r=0.4$","$r=0.3$"],loc="center",fontsize=16,frameon=False)
	gca().set_aspect("equal")
	axis("off")
