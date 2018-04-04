import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda
import pycuda.autoinit
import numpy
from pycuda.curandom import rand as curand



from pycuda.elementwise import ElementwiseKernel

lin_fun = ElementwiseKernel("float w2, float x1, float y1, float x2, float y2, float *I, float *J, float *t", 

    """
    t[i] = ((I[i]-x1)*(x2-x1)+(J[i]-y1)*(y2-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    
    if(t[i]<0.0){t[i]=0.0;}
    if(t[i]>1.0){t[i]=1.0;} 

    float xx; 
    float yy; 

    xx = x1+(x2-x1)*t[i]; 
    yy = y1+(y2-y1)*t[i]; 
    
    t[i] = (I[i]-xx)*(I[i]-xx)+(J[i]-yy)*(J[i]-yy); 
    
    if(t[i]<= 4*w2){
    t[i]=exp(-t[i]/(2.*w2));} 
    else{
    t[i]=0.0;}
    
    """)

"""
a = ones((100,100))
b = ones_like(a)
c = zeros_like(a)

a_gpu = gpuarray.to_gpu(a.astype(float32))
b_gpu = gpuarray.to_gpu(b.astype(float32))
c_gpu = gpuarray.to_gpu(c.astype(float32))

I,J = indices(shape(a))

I_gpu = gpuarray.to_gpu(I.astype(float32))
J_gpu = gpuarray.to_gpu(J.astype(float32))

w2 = 4.
x1= 50.
y1 = 50.
x2 = 75.
y2 = 75.

lin_fun(w2,x1,y1,x2,y2,I_gpu,J_gpu,c_gpu)

c = c_gpu.get()

imshow(c)
"""

