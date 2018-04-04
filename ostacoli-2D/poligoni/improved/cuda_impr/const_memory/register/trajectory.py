import matplotlib.colorbar as mcolorbar
D=1;tau=0.1
traettoria=[]
V=[]
for i in arange(0,9999):
	x,y=transpose(loadtxt("./posizione/dati_"+str(i)))
	traettoria=append(traettoria,[x[0],y[0]])
	vx,vy=transpose(loadtxt("./velocity/velocity_"+str(i)))
	V=append(V,[vx[0],vy[0]])
x=traettoria[0::2]
y=traettoria[1::2]
vx=V[0::2]
vy=V[1::2]
modulo=sqrt(vx**2+vy**2)
def stampa_traettoria():
	cond=slice(835,920,1)
	cond1=slice(835+1,920+1,1)
	#plot(x[cond],y[cond],"-",ms=6,mec="c",mew=2,c="g",lw=3)
	modulo=sqrt(vx**2+vy**2)[cond]
	m=min(modulo)
	M=max(modulo)
	norm=Normalize(m,M)
	norm.autoscale(modulo)
	f,ax1=subplots(1,1,figsize=(4,3),sharex='all',sharey='all')
	quiver(x[cond],y[cond],vx[cond1],vy[cond1],color=cm.jet(norm(modulo)),scale=100)
	fill_between([-.3,0],[1,1],[1.1,1.1],color="gray",alpha=0.5)
	gca().set_aspect("equal")
	xlim(-0.14,-0.012)
	ylim(0.93,1.005)
	xlabel("$x$")
	ylabel("$y$")
	tight_layout()
	xlabel("$x$",fontsize=16)
	ylabel("$y$",fontsize=16)

	xticks(np.arange(-0.14,-0.012, .04))
	yticks(np.arange(.93,1.005, .02))

	f.subplots_adjust(right=0.75)
	ax1 = f.add_axes([0.8, 0.26, 0.05, 0.48])
	mcolorbar.ColorbarBase(ax1, cmap='jet',norm=norm,orientation='vertical')

def crea_frames(start=1500,stop=3200):
	x_sh=x[start:stop:6]
	y_sh=y[start:stop:6]
	vx_sh=vx[start:stop:6]
	vy_sh=vy[start:stop:6]
	modulo_sh=modulo[start:stop:6]/sqrt(D/tau)
	m=min(modulo_sh)
	M=max(modulo_sh)
	norm=Normalize(m,M)
	norm.autoscale(modulo_sh)
	t=len(x_sh)
	f=figure(figsize=(4,3))
	xlim(-1.1,1.1)
	ylim(-1.1,1.1)
	xlabel("$x$",fontsize=10)
	ylabel("$y$",fontsize=10)
	tight_layout()
	xlabel("$x$",fontsize=16)
	ylabel("$y$",fontsize=16)
	gca().set_aspect("equal")
	ax=gca()
	fill_between([-1.1,1.1],-1.1,-1.,color="gray",alpha=0.4)
	fill_between([-1.1,1.1],1.,1.1,color="gray",alpha=0.4)
	fill_between([-1.1,-1.],-1.,1.,color="gray",alpha=0.4)
	fill_between([1.,1.1],-1.,1.,color="gray",alpha=0.4)

	f.subplots_adjust(right=0.75)
	ax1 = f.add_axes([0.8, 0.18, 0.05, 0.77])
	mcolorbar.ColorbarBase(ax1, cmap='jet',norm=norm,orientation='vertical')
	for i in arange(0,t):
		#plot(x[i],y[i],"o",ms=6,mec="c",mew=2,c="g",lw=3)
		ax.quiver(x_sh[i],y_sh[i],vx_sh[i],vy_sh[i],color=cm.jet(norm(modulo_sh[i])),scale=150)
		savefig("./frames/"+str(i))
		#wait= input("press something To continue")
		#print("ok andiamo avanti")
	#ffmpeg-f image2 -i ./frames/%d.png -sameq -r 25 movie.avi
	

