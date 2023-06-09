# self use
import numpy as np
from scipy.optimize import leastsq
from scipy import integrate
import matplotlib.pyplot as plt
T = 100
x = 1
N = 8
kb = 1.3806 * 10**-23
gamma = 1.19
v0 = 3.566*10**-10
K = 450 * 10**9
uu=1.66*10**-27
M=1
C124 = "J:\diamond\data.txt"
C123 = "J:\diamond\C-12-13-a"

x_data=[]
y_data=[]

with open(C124) as fdata:
    for line in fdata:
        data = line.strip().split('\t')
        x_data.append(float(data[0]))
        y_data.append((float(data[2])*0.5*0.529177249)**3)

V = y_data[0]*10**-30

def debyeintef(x):
    if x == 0.0:
        f = 0.0
    else:
        f = x**4 * np.exp(x)/(np.exp(x)-1)**2
    return f

def debye_cv(t,theta_D):
    cv, error = integrate.quad(debyeintef, 0, theta_D/t)
    return 9 * N * kb * (t/theta_D)**3 * cv

def alpha(t,para):
    D1,g1=para
    D=D1*M**0.02464
    g=g1*M**-0.2134
    alpha = (g * debye_cv(t,abs(D)) )/ K / V
    return alpha

def Volume(t,para):

    Volume, error = integrate.quad(alpha, 100, t,args=para)
    Volume = V * np.exp(Volume)*10**30
    return Volume

def error (para,t,y):
    result = []
    for t_val, y_val in zip(t,y):
        #print(t_val,y_val)
        result.append(Volume(t_val,para)-y_val)
    return result







para = leastsq(error,[1000,0.5],args=(x_data,y_data))
print(para[0])


y_fit = []
for xx in x_data:
    y_fit.append(Volume(xx,para[0]))

plt.figure
plt.plot(x_data,y_data,'r', label = "Ori")
plt.plot(x_data,y_fit,'-b',label ='fit')
plt.legend
plt.show()
