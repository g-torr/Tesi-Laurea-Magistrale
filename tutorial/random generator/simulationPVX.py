execfile("pycu_import.py")

sim_mod = SourceModule("""

#include <cuda.h>
#include <math.h>

#define PI 3.141592653589793f
#define SQ2 1.4142135623730951f

#define rc 2.5f

#define rc2 powf(rc,2)

#define Lx 20.0f

__device__ float rng_uni(unsigned int *state)
{
        // generates uniform ran num, 
        // keeps state val pointed by *state updated 

        unsigned int x = *state;

        x = x ^ (x >> 13);
        x = x ^ (x << 17);
        x = x ^ (x >> 5);

        *state = x;

        return (float) x / 4294967296;
}


__device__ void bm_trans(float& u1, float& u2)
{
        // transforms two given unif ran 
        // nums into two gauss ran nums

        float r = sqrtf(-2.0f * logf(u1));
        float phi = 2.0f * PI * u2;
        u1 = r * cosf(phi);
        u2 = r * sinf(phi);
}

__device__ void rng_gauss(float& u1, float& u2, unsigned int *state)
{
        // given two floats u1 and u2 it writes on them
        // two unif rand nums between 0 and 1 and transform     
        // them in gauss ran nums WARNING: variance is 1 !

        u1 = rng_uni(state);
        u2 = rng_uni(state);

        bm_trans(u1,u2);
}

__device__ float fx( float x)
{
    float myfx;

    //myfx = -x;

    myfx = x-x*x*x;

    return myfx;

}

__device__ float pbcx( float x)
{
    float xp = x;

    if( x < -Lx*0.5f)
    {
        xp = x+Lx;
    }

    if( x > Lx*0.5f)
    {
        xp = x-Lx;
    }

    return xp;
}

__device__ void RK2( float *x, float *y, float *nx, float *ny, float *Fx, float *Fy, unsigned int *state, float dt, float dt_sq, float D_sq, float tau)
{
    float x0 = *x;
    float y0 = *y;

    float w1 = 0.0f;
    float w2 = 0.0f;

    float nxi = *nx;
    float nyi = *ny;

    rng_gauss(w1,w2,state);

    w1 = SQ2*w1;
    w2 = SQ2*w2;

    *nx = nxi - (1.0f/tau)*nxi*dt + (D_sq/tau)*dt_sq*w1;
    *ny = nyi - (1.0f/tau)*nyi*dt + (D_sq/tau)*dt_sq*w2;

    *Fx = fx(x0);
    *Fy = fx(y0);

    *x =  x0  + dt*(*Fx) + 1.0*dt*(nxi) ;
    *y =  y0  + dt*(*Fy) + 1.0*dt*(nyi) ;


}


__global__ void dyn_kernel( float *x, float *y, float *nx, float *ny, float *Fx, float *Fy, unsigned int *states, float dt, float dt_sq, float D_sq, float tau, int N)
{
        // KERNEL: launches parallel stochastic evlution on threads

        int tid = threadIdx.x+blockIdx.x*blockDim.x;
        // WARNING: better to be sobstituted with absolute thread coord...

        if (tid < N)
        {
                // float u1 = 0.;
                // float u2 = 0.;

                RK2( &x[tid], &y[tid], &nx[tid], &ny[tid], &Fx[tid], &Fy[tid], &states[tid], dt, dt_sq, D_sq, tau);
        }
}

""")


import os


class Simulation:

    def __init__( self, NB, NT, D , tau ):
        #---------------------------------
        self.NB = int(NB)
        self.NT = int(NT)
        self.N = int32(NB*NT)
        #----------------------------------
        self.D = float32(D)
        self.tau = float32(tau)
        #----------------------------------
        self.D_sq = float32( sqrt(D))
        #----------------------------------
        self.initialize_dtoh()
        #----------------------------------
        self.evolve = sim_mod.get_function("dyn_kernel")
        print '\n... simulation object initialized!'
        self.print_info()


    def simulate( self, dt, Ttot, Teq, Tsave):
        #-----------------------------
        self.dt = float32(dt)
        self.dt_sq = float32(sqrt(dt))
        self.Ttot = int(Ttot)
        self.Teq = int(Teq)
        self.Tsave = int(Tsave)
        self.data = []
        #------------------------------
        self.print_dyn_info()
        #------------------------------
        self.hi = zeros(99)
        vv = 5.*(self.D/self.tau)
        bns = linspace(-vv,vv,100)
        self.F = []
        self.F2 = []
        self.nsaves = 0
        self.xu = []
        self.ig = []
        self.xu2 = []
        self.ig2 = []
        self.fxe = []
        self.Q = []
        self.S = []
        tup = int(float(self.Ttot)/10000.)
        if tup==0:
            tup = 1
        self.pvx = None
        print 'update time  = '+str(tup)+'\n'
        
        for t in range(0,self.Ttot):
            
            prog = int(100.*t/float(self.Ttot))
            if (t%tup )==0:
                print 'simulation progress =  '+str(prog)+' % \r',

            self.evolve( self.x_gpu, self.y_gpu, self.nx_gpu, self.ny_gpu, \
                    self.fx_gpu, self.fy_gpu, self.states_gpu, \
                    self.dt, self.dt_sq, self.D_sq, self.tau, self.N, \
                    block=(self.NT, 1, 1), grid=(self.NB, 1) )

            if t > self.Teq:
                if self.pvx ==None:
                    self.get_all_data_dtoh()
                    xstd = sqrt( (var(self.x)+var(self.y))/2. )
                    vstd = sqrt( (var(self.nx+self.fx)+var(self.ny+self.fy))/2. )
                    print "\n setting hist with stdx = "+str(xstd)+" stdv = "+str(vstd)+"\n"
                    xvbns = [linspace(-5*xstd,5*xstd,200),linspace(-5*vstd,5*vstd,202)]                                      
                    hx,bxx,bvv = histogram2d(self.x,self.nx+self.fx,bins = xvbns)
                    """figure()
                    hist(self.x)
                    figure()
                    hist(self.nx + self.fx)
                    figure()
                    imshow(hx)"""
                    self.pvx = zeros_like(hx)
                if (t % self.Tsave == 0):
                    self.nsaves += 1
                    self.get_all_data_dtoh()
                    h,b = histogram( self.x, bins = bns) 
                    self.hi += h
                    self.fxe += [(mean(self.x*self.fx)+mean(self.y*self.fy))/2. ]
                    self.xu += [(mean(self.x*self.nx)+mean(self.y*self.ny))/2. ]
                    self.ig += [(mean(self.D/(1.+self.tau))+mean(self.D/(1.+self.tau)) )/2.]
                    self.xu2 += [(var(self.x*self.nx)+var(self.y*self.ny))/2. ]
                    self.ig2 += [(var(self.D/(1.+self.tau))+var(self.D/(1.+self.tau)) )/2.]
                    h,b = histogram( self.y, bins = bns)
                    self.hi += h
                    self.F += [mean(abs(self.fx))]
                    self.F += [mean(abs(self.fy))]
                    self.F2 += [mean(self.fx**2)]
                    self.F2 += [mean(self.fy**2)]
                    hx,bxx,bvv = histogram2d(self.x,self.nx+self.fx,bins = xvbns)                    
                    hy,bxx,bvv = histogram2d(self.y,self.ny+self.fy,bins = xvbns)
                    self.pvx += (hx+hy)/2.
                    #--------------------------------------------------------------------------------------------------------------
                    QX = mean( self.Gamma(self.x)*self.theta(self.x)/self.tau -  self.Gamma(self.x)*(self.fx+self.nx)**2/self.tau )
                    #QX = mean( self.comq(self.x,self.fx+self.nx) )
                    #
                    QY = mean( self.Gamma(self.y)*self.theta(self.y)/self.tau -  self.Gamma(self.y)*(self.fy+self.ny)**2/self.tau )         
                    #QY = mean( self.comq(self.y,self.fy+self.ny) )
                    #
                    #SX = mean( ( self.Gamma(self.x)*self.theta(self.x)/self.tau -  self.Gamma(self.x)*(self.fx+self.nx)**2/self.tau )\
                            #        /self.theta(self.x) )                             
                    SX = mean( self.comqt(self.x,self.fx+self.nx) )
                    #
                    #SY = mean( ( self.Gamma(self.y)*self.theta(self.y)/self.tau -  self.Gamma(self.y)*(self.fy+self.ny)**2/self.tau )\
                            #/self.theta(self.y))                    
                    SY = mean( self.comqt(self.y,self.fy+self.ny) )
                    self.Q += [(QX+QY)/2.]
                    self.S += [(SX+SY)/2.]
        figure()
        subplot(2,1,1)
        plot(self.Q)
        semilogx()
        ylabel("Q")
        subplot(2,1,2)
        plot(self.S)
        semilogx()
        ylabel("S")
        self.dQ = array(self.Q).std()        
        self.dS = array(self.S).std()
        self.Q = array(self.Q).mean()
        self.S = array(self.S).mean()
        print "#############################################"
        print "total Q = "+str(self.Q)+" +/- "+str(self.dQ)+" +/- "+str(self.dQ/sqrt(self.nsaves))
        print "total S = "+str(self.S)+" +/- "+str(self.dS)+" +/- "+str(self.dS/sqrt(self.nsaves))
        print "#############################################"
        self.pvx = self.pvx/self.nsaves
        self.bxx = bxx
        self.bvv = bvv
        self.bx = (b[1:]+b[:-1])/2.
        self.hi = self.hi/(2.*self.nsaves)
        print "\n "
        self.xu = mean(array(self.xu))
        self.ig = mean(array(self.ig))
        self.fxe = mean(array(self.fxe))
        #print "mean ig = "+str(mean(array(self.ig)))+" +/- "+str(sqrt(mean(array(self.ig2))))
        #figure()
        #plot( self.bx, self.hi, '-o')
        self.histo = transpose(array([self.bx, self.hi]))
        self.F = mean( self.F)
        self.F2 = mean( self.F2)
        print 'num of samples : '+str(self.nsaves)
        print '\n ... simulation finished!'

    def comq(self,x,v):
        D = self.D
        tau = self.tau
        ph2 = self.phi2(x)
        return ( D-(v**2)*tau*(1.+tau*ph2 ) )/tau**2

    def comqt(self,x,v):
        D = self.D
        tau = self.tau
        ph2 = self.phi2(x)
        return (1+tau*ph2)*(D -(v**2)*tau - ph2*(v*tau)**2)/(D*tau)
    

    def phi2(self,x):
        return -1+3*x**2

    def Gamma(self,x):
        return 1.-self.tau*(1.-3*x**2)
    
    def theta(self,x):
        Tb = self.D/self.tau
        return Tb/self.Gamma(x)

    def plteo(self):
        self.hi = self.hi/trapz(self.hi,self.bx)
        figure()
        plot(self.bx,self.hi,'-',lw=2,label = "Numerical")
        x = linspace(self.bx.min(),self.bx.max(),1000)
        u = (x**4)/4. - (x**2)/2.
        u1 = (x**3) - x
        u2 = 3*(x**2) -1.
        p0 = exp(-u/self.D) 
        p1 = exp(-self.tau*(u1**2)/(2.*self.D))
        p2 = abs(1.+ self.tau *u2)
        pu = p0*p1*p2
        pu = pu/trapz(pu,x)
        pumx = max(pu)
        plot(x,pu,'-m',label = "UCNA")
        xm = -1./sqrt(3.)
        xp = 1./sqrt(3.)
        plot([xm,xm],[1e-3,1.1*pumx],'--k')
        plot([xp,xp],[1e-3,1.1*pumx],'--k')
        ylim(0,1.1*pumx)
        p2 = where(u2>0,1.+ self.tau *u2,1./(1-self.tau*u2))
        pu = p0*p1*p2
        pu = pu/trapz(pu,x)
        plot(x,pu,'-r',label="Hybrid")
        p01 = exp( cumsum(-u1/(1-self.tau*u2))*diff(x)[0] )
        p2 = 1./abs(1.- self.tau *u2)
        pu = p01*p2
        pu = pu/trapz(pu,x)
        plot(x,pu,'-g',label = "Small $\\tau$")
        pb = p0/trapz(p0,x)
        #plot(x,pb,'-k',label = "Boltzmann")
        legend()

    def get_all_data_dtoh( self):
        cuda.memcpy_dtoh( self.x, self.x_gpu) 
        cuda.memcpy_dtoh( self.y, self.y_gpu)
        cuda.memcpy_dtoh( self.nx, self.nx_gpu) 
        cuda.memcpy_dtoh( self.ny, self.ny_gpu)
        cuda.memcpy_dtoh( self.fx, self.fx_gpu)
        cuda.memcpy_dtoh( self.fy, self.fy_gpu)
        cuda.memcpy_dtoh( self.states, self.states_gpu)

    def save_data( self, folnam):
        
        self.folnam = folnam
        
        os.makedirs(self.folnam)

        self.simpar = array([self.NB, self.NT, self.N, self.D, self.tau, \
                self.dt, self.Ttot, self.Teq, self.Tsave, self.nsaves]) 

        np.save( self.folnam+'/simpar.npy', self.simpar)
        np.save( self.folnam+'/histo.npy', self.histo)
        np.save( self.folnam+'/F_F2.npy', array([self.F, self.F2]))



    def print_info(self):
        print '\n -----INFOS-----------'
        print 'NB = '+str(self.NB)
        print 'NT = '+str(self.NT)
        print 'N = '+str(self.N)
        print '----------------------'
        print 'D = '+str(self.D)
        print 'tau = '+str(self.tau)
        print '------------------------\n'
        
    
    def print_dyn_info(self):
        print '\n -----INFOS-----------'
        print 'dt = '+str(self.dt)
        print 'Ttot = '+str(self.Ttot)+' = '+str(self.Ttot*self.dt)
        print 'Teq = '+str(self.Teq)+' = '+str(self.Teq*self.dt)
        print 'Tsave = '+str(self.Tsave)+' = '+str(self.Tsave*self.dt)
        print '------------------------\n'

    
    def initialize_dtoh(self):

        N = self.N
        D = self.D
        tau = self.tau

        #-----------------------------------
        states = arange(1,2**32, int(2.**32/N)).astype(uint32)
        #states = ones(N).astype(uint32)
        states_gpu = cuda.mem_alloc(states.nbytes)
        cuda.memcpy_htod(states_gpu, states)
        self.states = states
        self.states_gpu = states_gpu
        
        #------------------------------------
        nx = (sqrt(D/tau)*randn(N)).astype(float32)
        nx_gpu = cuda.mem_alloc(nx.nbytes)
        cuda.memcpy_htod(nx_gpu, nx)
        self.nx = nx
        self.nx_gpu = nx_gpu
        
        ny = (sqrt(D/tau)*randn(N)).astype(float32)
        ny_gpu = cuda.mem_alloc(ny.nbytes)
        cuda.memcpy_htod(ny_gpu, ny)
        self.ny = ny
        self.ny_gpu = ny_gpu
        
        #-----------------------------------
        x = (2.5* ones(N)).astype(float32)
        x_gpu = cuda.mem_alloc(x.nbytes)
        cuda.memcpy_htod(x_gpu, x)
        self.x = x
        self.x_gpu = x_gpu
        
        y = (-2.5*ones(N)).astype(float32)
        y_gpu = cuda.mem_alloc(y.nbytes)
        cuda.memcpy_htod(y_gpu, y)
        self.y = y
        self.y_gpu = y_gpu
        
        #--------------------------------
        fx = zeros(N).astype(float32)
        fx_gpu = cuda.mem_alloc(fx.nbytes)
        cuda.memcpy_htod(fx_gpu, fx) 
        self.fx = fx
        self.fx_gpu = fx_gpu
        
        fy = zeros(N).astype(float32)
        fy_gpu = cuda.mem_alloc(fy.nbytes)
        cuda.memcpy_htod(fy_gpu, fy)
        self.fy = fy
        self.fy_gpu = fy_gpu


    
   




def scan_param_space():
    NB = 16
    NT = 128
    dt = 1e-4
    Tsave = 1e3
    Teq = 5e6
    Ttot = 6e6
    for D in linspace(0.4,100.,5):
        for tau in array([0.1,0.325,1.]):
            folnam = 'data_'+str(D)+'_'+str(tau)
            so = Simulation(  NB, NT, D , tau )
            so.simulate( dt, Ttot, Teq, Tsave)
            so.save_data( folnam)

def scan_param_space_new():
    NB = 16
    NT = 128
    dt = 1e-4
    Tsave = 1e3
    Teq = 5e6
    Ttot = 6e6
    #--------------------------------------
    D = 0.4
    tau = 0.1
    folnam = 'new_data_'+str(D)+'_'+str(tau)
    so = Simulation(  NB, NT, D , tau )
    so.simulate( dt, Ttot, Teq, Tsave)
    so.save_data( folnam)
    #--------------------------------------
    D = 100.0
    tau = 1.0
    folnam = 'new_data_'+str(D)+'_'+str(tau)
    so = Simulation(  NB, NT, D , tau )
    so.simulate( dt, Ttot, Teq, Tsave)
    so.save_data( folnam)


def control_noise():
    NB = 4
    NT = 512
    dt = 1e-4
    #Tsave = 1e3
    #Teq = 5e6
    #Ttot = 6e6
    figure()
    for D in linspace(0.4,100.,2):
        for tau in array([0.1,0.325,1.]):
            Tsave = tau/1e-4
            Teq = 10.*tau/1e-4
            Ttot = 100.*tau/1e-4
            #folnam = 'data_'+str(D)+'_'+str(tau)
            so = Simulation(  NB, NT, D , tau )
            nx0 = so.nx
            so.simulate( dt, Ttot, Teq, Tsave)
            nx1 = so.nx
            print 'checking difference in noise = '+str(nx1-nx0)
            #so.save_data( folnam)
            h, b = histogram(so.nx,bins=100,range=[-3*sqrt(D/tau),3*sqrt(D/tau)])
            bx = (b[1:]+b[:-1])/2.
            plot(bx,h)
            title('GCpin')





