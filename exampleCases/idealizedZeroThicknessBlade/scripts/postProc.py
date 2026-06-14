#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 21:13:12 2025

@author: adekola
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

filePath  = '/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/cases/'

Nt = 360
Nr = 3
theta = np.linspace(0, 2*np.pi, Nt)
def massAvgQty(Pt, Uz, Nt, Nr, stagger_deg=10, axialBladeIn=0, axialBladeOut=1):
    
    # Rest of your code remains the same
    radius = np.sqrt(Pt[:,0]**2 + Pt[:,1]**2)
    hub = min(radius)*1.001
    cas = max(radius)*0.999
    spanFrac = np.linspace(0,1,Nr)
    radi = hub + spanFrac*(cas - hub)
    theta = np.linspace(0, 2*np.pi, Nt) 
    index = 0
    data = np.zeros([Nr*Nt, 2])
    X = np.zeros((Nt, Nr))
    Y = np.zeros((Nt, Nr))
    for b in range(Nr):
        for c in range(Nt):   
            thetaOffset = (axialBladeOut - axialBladeIn) * np.tan(np.deg2rad(stagger_deg))
            # Convert arc length to angle at this radius
            angle = thetaOffset / radi[b]
            # Total angle at outlet = base theta + stagger rotation
            theta_total = theta[c] + angle
            
            X[c, b] = radi[b] * np.cos(theta_total)
            Y[c, b] = radi[b] * np.sin(theta_total)
            # data[index] = (radi[b] * np.sin(theta[c])),(radi[b] * np.cos(theta[c]))
            # X[c,b] = (radi[b] * np.sin(theta[c]))
            # Y[c,b] = (radi[b] * np.cos(theta[c])) 
            index +=1
    Area = np.zeros(Nr*Nt) #A=0.5*​(router2​−rinner2​)Δθ
    dr = np.diff(radi)[0]  # Radial spacing
    dtheta = np.diff(theta)[0]  # Angular spacing
    count = 0
    for c in range(Nr):
        for d in range(Nt):
            r = radi[c]
            # Define cell boundaries
            rInner = r - 0.5 * dr
            rOuter = r + 0.5 * dr
            # Handle boundary cases
            if c == 0:  # Inner boundary
                rInner = hub
            elif c == Nr - 1:  # Outer boundary  
                rOuter = cas       
            # Ensure positive radii
            rInner = max(rInner, 0)
            rOuter = max(rOuter, rInner + 1e-10)      
            # Correct annular sector area
            Area[count] = 0.5 * (rOuter**2 - rInner**2) * dtheta
            count += 1
    print(Area.shape)
    newPt = griddata((Pt[:,0], Pt[:,1]), Pt[:,3], (X, Y), method='cubic')
    newUz = griddata((Uz[:,0], Uz[:,1]), Uz[:,5], (X, Y), method='cubic')
    massAvgPt = np.zeros(Nr)
    for i in range(Nr):
        AvgPt = np.zeros(Nt)
        mDot= np.zeros(Nt)
        for j in range(Nt):
            # print(f"j={j}, i={i}")
            # print(f"newPt[j,i] shape: {np.shape(newPt[j,i])}, value: {newPt[j,i]}")
            # print(f"Area[j+Nt*i] shape: {np.shape(Area[j+Nt*i])}, value: {Area[j+Nt*i]}")
            # print(f"newUz[j,i] shape: {np.shape(newUz[j,i])}, value: {newUz[j,i]}")
            # print(f"newRho[j,i] shape: {np.shape(newRho[j,i])}, value: {newRho[j,i]}")
            AvgPt[j] = newPt[j,i]*Area[j+Nt*i]*newUz[j,i]
            mDot[j] = Area[j+Nt*i]*newUz[j,i]
        massAvgPt[i] = sum(AvgPt) / sum(mDot)
    return massAvgPt, spanFrac, newPt, X, Y

rCase1Pt = np.loadtxt(filePath + '/baseLine/postProcessing/rotorOutlet/1027/totalP_zNormal.raw')
rCase1AlphaFlow = np.loadtxt(filePath + '/baseLine/postProcessing/rotorOutlet/1027/alphaFlow_zNormal.raw')
rCase1Uz = np.loadtxt(filePath + '/baseLine/postProcessing/rotorOutlet/1027/U_zNormal.raw')

rAvgCasePt1, rSpanFrac, rCase1Ptnew, X, Y = massAvgQty(rCase1Pt, rCase1Uz, Nt, Nr)
rAvgCase1, _, rCase1AlphaFlownew,_,_ = massAvgQty(rCase1AlphaFlow, rCase1Uz, Nt, Nr)

plt.plot(theta, rCase1AlphaFlownew[:,1], 'k')
plt.axis('equal')
#%%
# plt.plot(rAvgCase1,rSpanFrac, 'k')
plt.plot(rCase1Pt[:,0], rCase1Pt[:,1], 'k.')
plt.plot(X,Y,'r.')
plt.axis('equal')









