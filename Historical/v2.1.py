#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
#from qutip import *
import matplotlib.pyplot as plt
from scipy.linalg import expm


# In[2]:


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


# In[3]:


class CliffordGate:
    def __init__(self):
        self.angles = [ (1.3e-5, 0, 0), 
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
                        (np.pi*2/3, np.arccos(1/np.sqrt(3)), np.pi*5/4)   ]
    
    def getgate_angle(self, gamma, theta, phi):
        self.target = np.array([[np.cos(gamma/2)+(0+1j)*np.sin(gamma/2)*np.cos(theta), 0, (0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0-1j)*phi)],
                                [0, 1, 0],
                                [(0-1j)*np.sin(gamma/2)*np.sin(theta)*np.exp((0+1j)*phi), 0, np.cos(gamma/2)-(0+1j)*np.sin(gamma/2)*np.cos(theta)]])  
        
    def getgate(self, n):
        self.gamma = self.angles[n-1][0]
        self.theta = self.angles[n-1][1]
        self.phi = self.angles[n-1][2]
        self.target = np.array([[np.cos(self.gamma/2)+(0+1j)*np.sin(self.gamma/2)*np.cos(self.theta), 0, (0-1j)*np.sin(self.gamma/2)*np.sin(self.theta)*np.exp((0-1j)*self.phi)],
                                [0, 1, 0],
                                [(0-1j)*np.sin(self.gamma/2)*np.sin(self.theta)*np.exp((0+1j)*self.phi), 0, np.cos(self.gamma/2)-(0+1j)*np.sin(self.gamma/2)*np.cos(self.theta)]])  
    
    def randtarget(self):
        # Generalize a random integer in the range of 1 to 24
        randnum = np.random.randint(low=1, high=24) # high will be 24 in the later version
        # Pick out the corresponding value of angles
        self.gamma = self.angles[randnum-1][0]
        self.theta = self.angles[randnum-1][1]
        self.phi = self.angles[randnum-1][2]
        Miu = (2*np.pi-2*self.gamma)/(2*np.pi)
        self.time = 2*np.pi/v1*np.sqrt(1-Miu**2)+10e-6
        self.eta=-2*np.pi*Miu/self.time
        # Utarget1=[cos(ThetaA1) 0 sin(ThetaA1)*exp(-1i*Fei1);0 1 0;sin(ThetaA1)*exp(1i*Fei1) 0 -cos(ThetaA1)];
        # Factor  np.exp((0-1j)*(self.gamma/2))* omitted
        self.target = np.array([[np.cos(self.gamma/2)+(0+1j)*np.sin(self.gamma/2)*np.cos(self.theta), 0, (0-1j)*np.sin(self.gamma/2)*np.sin(self.theta)*np.exp((0-1j)*self.phi)],
                                [0, 1, 0],
                                [(0-1j)*np.sin(self.gamma/2)*np.sin(self.theta)*np.exp((0+1j)*self.phi), 0, np.cos(self.gamma/2)-(0+1j)*np.sin(self.gamma/2)*np.cos(self.theta)]])  


# In[4]:


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
    if gammaq == 0:
        thetaq =0
    else:
        thetaq = np.real(np.arccos(((temp1-temp9)/2)/((0+1j)*np.sin(gammaq/2))))
    if thetaq==0 or gammaq==0:
        phiiq = 0
    else:
        phiiq = -np.arcsin(np.real(temp3)/np.sin(gammaq/2)/np.sin(thetaq))
        if np.allclose(UI(gammaq,thetaq,phiiq), U) == False:
            phiiq = np.pi - phiiq
    Miu = (2*np.pi-2*gammaq)/(2*np.pi)
    time = 2*np.pi/v1*np.sqrt(1-Miu**2)
    eta = -2*np.pi*Miu/time
    return gammaq, thetaq, phiiq, time, eta


# In[5]:


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
v1=80.7954


# In[6]:


# Initial states
intialstate = zerostate
#intialstate =  np.multiply(1/np.sqrt(2), np.add(zerostate, onestate))
initialdm = np.matmul(intialstate, np.matrix(intialstate).getH())
dm1 = initialdm
finalatomstate1 = zerostate
fidm1 = np.matmul(finalatomstate1, np.matrix(finalatomstate1).H)
finalatomstate2 = estate
fidm2 = np.matmul(finalatomstate2, np.matrix(finalatomstate2).H)
finalatomstate3 = onestate
fidm3 = np.matmul(finalatomstate3, np.matrix(finalatomstate3).H)

# Max Rabi frequency
# Dots drawn per period每一个周期画的点数
dotsPerPeriod = 1000;
numberOfPeriod = 1;
# Define time and setp
#U = np.identity(3)
# Define time and setp
#timestep = tao/dotsPerPeriod;  # Draw dotsPerPeriod dots per period
#time = (m+1)*numberOfPeriod*dotsPerPeriod;  # Draw numberOfPeriod periods


# In[7]:


m = 10
Clifford = CliffordGate()
U_final = np.identity(3)


# In[8]:


for k in range(1, m+2, 1):
    if k == m+1:
        # Deal with the last gate - the recover gate
        U_recover = np.matrix(U_final).H
        gamma, theta, phi, time, eta = get_para(U_recover)
        #print(k, ":", gamma, theta, phi, time, eta, "\n")
    else:   
        Clifford.randtarget()
        U = Clifford.target
        gamma = Clifford.gamma
        theta = Clifford.theta
        phi = Clifford.phi
        eta = Clifford.eta
        time = Clifford.time
        U_final = np.matmul(U, U_final)
        #print(k, ":", gamma, theta, phi, time, eta, "\n", U_final, "\n")
    timestep = time/dotsPerPeriod
    for n in range(1, dotsPerPeriod+1, 1):        
        t = (n-1)*timestep   
        PHI = eta*t
        Omega = v1*np.exp((0-1j)*PHI)
        OmegaP = Omega*np.sin(theta/2)
        OmegaS = Omega*np.cos(theta/2)*np.exp((0+1j)*phi)
        Hamiltonian1 = (OmegaP*A0e+np.conj(OmegaP)*Ae0+OmegaS*A1e+np.conj(OmegaS)*Ae1)/2
        dt = timestep
        #U = np.matmul(expm((0-1j)*Hamiltonian1*dt), U)
        Psi = np.matmul(expm((0-1j)*Hamiltonian1*dt), Psi)
        k1 = solve_master(Hamiltonian1,dm1,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
        k2 = solve_master(Hamiltonian1,dm1+0.5*timestep*k1,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
        k3 = solve_master(Hamiltonian1,dm1+0.5*timestep*k2,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
        k4 = solve_master(Hamiltonian1,dm1+timestep*k3,Gamma0,Gamma1,Gamma2,kappa0,kappa1,kappa2)
        dm1 = dm1+(timestep/6)*(k1+2*k2+2*k3+k4)
        dm1 = 0.5*(dm1+np.matrix(dm1).H)
        dm1 = dm1/(np.trace(dm1))
        rho = dm1
        
fidelity = np.real(np.trace(np.matmul(rho, fidm1)))


# In[9]:


fidelity


# In[ ]:




