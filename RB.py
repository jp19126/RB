#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''
-0.701184236509162 + 0.709973051882869j,      -0.0325875156346839 - 0.0328255976389318j,    0.0325875156346838 - 0.0328255976389317j
0.000291958531064516 - 0.0462534409077854j,   0.00113618408451785 - 0.00807218387098761j,   0.998863815915481 - 0.00807218387098759j
0.000291958531064545 - 0.0462534409077857j,   -0.998863815915482 - 0.00807218387098761j,    0.00113618408451833 - 0.00807218387098760j
'''


# In[1]:


import numpy as np
#from qutip import *
import matplotlib.pyplot as plt
from scipy.linalg import expm
import time as ti


# In[2]:


Cliffords = [(1.3e-5, 0, 0), 
             (np.pi, np.pi/2, 0),
             (np.pi, np.pi/2, np.pi/2), 
             (np.pi, 0, 0),
             (np.pi, np.pi*3/4, 0),
             (np.pi, np.pi/4, 0),
             (np.pi, np.pi/4, np.pi/2),
             (np.pi, np.pi*3/4, np.pi/2), 
             (np.pi, np.pi/2, np.pi/4),
             (np.pi, np.pi/2, np.pi*3/4),
             (np.pi/2, np.pi/2, 0),
             (np.pi/2, np.pi/2, np.pi),
             (np.pi/2, np.pi/2, np.pi/2),
             (np.pi/2, np.pi/2, np.pi*3/2),
             (np.pi/2, 0, 0),
             (np.pi/2, np.pi, 0),
             (np.pi*2/3, np.arccos(-1/np.sqrt(3)), np.pi/4),
             (np.pi*2/3, np.arccos(1/np.sqrt(3)), -np.pi/4),
             (np.pi*2/3, np.arccos(1/np.sqrt(3)), np.pi*3/4),
             (np.pi*2/3, np.arccos(-1/np.sqrt(3)), np.pi*5/4),
             (np.pi*2/3, np.arccos(1/np.sqrt(3)), np.pi/4),
             (np.pi*2/3, np.arccos(-1/np.sqrt(3)), np.pi*3/4),
             (np.pi*2/3, np.arccos(-1/np.sqrt(3)), -np.pi/4), 
             (np.pi*2/3, np.arccos(1/np.sqrt(3)), np.pi*5/4)]


# In[3]:


def solve_master(Hami,dm,decoop1,decoop2,decoop3,rate1,rate2,rate3):
    down1 = decoop1
    down2 = decoop2
    down3 = decoop3
    up1 = np.matrix(down1).H
    up2 = np.matrix(down2).H
    up3 = np.matrix(down3).H
    temp1 = (0-1j)*(np.dot(Hami, dm) - np.dot(dm, Hami))
    temp2 = 2*np.dot(np.dot(down1, dm), up1)-np.dot(np.dot(up1, down1), dm)-np.dot(np.dot(dm, up1), down1)
    temp3 = 2*np.dot(np.dot(down2, dm), up2)-np.dot(np.dot(up2, down2), dm)-np.dot(np.dot(dm, up2), down2)
    temp4 = 2*np.dot(np.dot(down3, dm), up3)-np.dot(np.dot(up3, down3), dm)-np.dot(np.dot(dm, up3), down3)
    out =  temp1+(np.multiply(rate1, temp2)+np.multiply(rate2, temp3)+np.multiply(rate3, temp4))/2
    return out


# In[4]:


def getgate(gamma, theta, phi):
    target = np.array([[np.cos(gamma/2)+(0+1j)*np.sin(gamma/2)*np.cos(theta), 0, (0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0-1j)*phi)],
                            [0, 1, 0],
                            [(0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0+1j)*phi), 0, np.cos(gamma/2)-(0+1j)*np.sin(gamma/2)*np.cos(theta)]])  
    return target


# In[5]:


def randClifford():
    randnum = np.random.randint(low=1, high=25) # Generalize a random integer in the range of 1 to 24 
    # Pick out the corresponding value of angles
    gamma = Cliffords[randnum-1][0]
    theta = Cliffords[randnum-1][1]
    phi = Cliffords[randnum-1][2]
    Miu = (2*np.pi-2*gamma)/(2*np.pi)
    time = 2*np.pi/v1*np.sqrt(1-Miu**2)+10e-6
    eta = -2*np.pi*Miu/time
    # Factor  np.exp((0-1j)*(gamma/2))* omitted
    target = np.array([[np.cos(gamma/2)+(0+1j)*np.sin(gamma/2)*np.cos(theta), 0, (0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0-1j)*phi)],
                       [0, 1, 0],
                       [(0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0+1j)*phi), 0, np.cos(gamma/2)-(0+1j)*np.sin(gamma/2)*np.cos(theta)]]) 
    
    time = 2*np.pi/v1*np.sqrt(1-Miu**2)+10e-6
    eta = -2*np.pi*Miu/time
    return target, gamma, theta, phi, time, eta


# In[6]:


def get_para(U):
    temp1 = U[0,0]
    temp9 = U[2,2]
    temp3 = U[0,2]
    temp7 = U[2,0]
    def UI(gamma, theta, phi):
        UII = np.array([[np.cos(gamma/2)+(0+1j)*np.sin(gamma/2)*np.cos(theta), 0, (0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0-1j)*phi)],
                       [0, 1, 0],
                       [(0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0+1j)*phi), 0, np.cos(gamma/2)-(0+1j)*np.sin(gamma/2)*np.cos(theta)]])  
        return UII
    gammaq = np.real(2*np.arccos((temp1+temp9)/2))
#     if np.allclose(gammaq, 2*np.pi):
#         gammaq = 0
    if gammaq == 0:
        thetaq = 0
    else:
        thetaq = np.real(np.arccos(((temp1-temp9)/2)/((0+1j)*np.sin(gammaq/2))))   
    if thetaq<1e-4 or gammaq<1e-4:
        phiiq = 0
    else:
        temp3_real = np.real(temp3)
        if np.allclose(temp3_real, 0):
            temp3_real = 0
        arg = temp3_real/np.sin(gammaq/2)/np.sin(thetaq)
        if arg>1 and np.allclose(arg, 1):
            arg = 1
        elif arg<-1 and np.allclose(arg, -1):
            arg = -1
        phiiq = -np.arcsin(arg)
        if np.allclose(UI(gammaq,thetaq,phiiq), U) == False:
            phiiq = np.pi - phiiq
    Miu = (2*np.pi-2*gammaq)/(2*np.pi)
    time = 2*np.pi/v1*np.sqrt(1-Miu**2)+10e-6
    eta = -2*np.pi*Miu/time
    return gammaq, thetaq, phiiq, time, eta


# In[7]:


identy = np.identity(3)
zerostate = [[1], [0], [0]]
estate = [[0], [1], [0]]
onestate = [[0], [0], [1]]
Psi = np.identity(3)
# Define the operator
A0e = np.dot(zerostate, np.transpose(estate))
Ae0 = np.dot(estate, np.transpose(zerostate))
Ae1 = np.dot(estate, np.transpose(onestate))
A1e = np.dot(onestate, np.transpose(estate))
A00 = np.dot(zerostate, np.transpose(zerostate))
Aee = np.dot(estate, np.transpose(estate))
A11 = np.dot(onestate, np.transpose(onestate))
# Define the decay operator
Gamma0 = A0e+A1e
Gamma1 = 2*Aee-A00-A11
Gamma2 = np.identity(3)
# Decay value
kappa0 = 2*np.pi*4/1000; # Cavity decay
kappa1 = 2*np.pi*4/1000; # Effective atom decay 16ns
kappa2 = 0;
# Area of pulse
v1 = 80.7954


# In[8]:


# Initial states
intialstate = zerostate
initialdm = np.matmul(intialstate, np.matrix(intialstate).getH())
finalatomstate1 = zerostate
fidm = np.matmul(finalatomstate1, np.matrix(finalatomstate1).H)

# Dots drawn per period每一个周期画的点数
dotsPerPeriod = 1000;
numberOfPeriod = 1;


# In[33]:


# Maximum number of Clifford gates
M = np.zeros(26, dtype = int)
for n in range(1, 27):
    M[n-1] = 4*(n-1)+1
# Number of iteration to get the average fidelity
K = 40


# In[10]:


fidelity = np.zeros(K)
fid_avg = np.zeros(len(M))


# In[24]:


start = ti.time()

for m in M:
    fidelity = np.zeros(K)
    for k in range (1, K+1, 1):
        dm1 = initialdm
        U_final = np.identity(3)
        for i in range(1, m+2, 1):
            if i == m+1:
                # Deal with the last gate - the recover gate
                U_recover = np.matrix(U_final).H
                gamma, theta, phi, time, eta = get_para(U_recover)
            else:   
                U, gamma, theta, phi, time, eta = randClifford()
                U_final = np.matmul(U, U_final)
            timestep = time/dotsPerPeriod
            for n in range(1, dotsPerPeriod+1, 1):        
                t = (n-1)*timestep   
                PHI = eta*t
                Omega = v1*np.exp((0-1j)*PHI)
                OmegaP = Omega*np.sin(theta/2)
                OmegaS = Omega*np.cos(theta/2)*np.exp((0+1j)*phi)
                Hamiltonian1 = (OmegaP*A0e+np.conj(OmegaP)*Ae0+OmegaS*A1e+np.conj(OmegaS)*Ae1)/2
                dt = timestep
                # U = np.matmul(expm((0-1j)*Hamiltonian1*dt), U)
                # Psi = np.matmul(expm((0-1j)*Hamiltonian1*dt), Psi)
                k1 = solve_master(Hamiltonian1,dm1,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
                k2 = solve_master(Hamiltonian1,dm1+0.5*timestep*k1,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
                k3 = solve_master(Hamiltonian1,dm1+0.5*timestep*k2,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
                k4 = solve_master(Hamiltonian1,dm1+timestep*k3,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
                dm1 = dm1+(timestep/6)*(k1+2*k2+2*k3+k4)
                dm1 = 0.5*(dm1+np.matrix(dm1).H)
                rho = dm1/(np.trace(dm1))
        fidelity[k-1] = np.real(np.trace(np.matmul(rho, fidm)))
    fid_avg[int((m-1)/4)] = np.mean(fidelity)
    
end = ti.time()


# In[ ]:


fid_avg


# In[ ]:




