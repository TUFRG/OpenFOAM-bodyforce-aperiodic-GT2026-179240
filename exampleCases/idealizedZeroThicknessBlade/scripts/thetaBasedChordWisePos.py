#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 21:39:06 2025

@author: adekola
"""

import numpy as np

def cart2pol(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.atan2(y,x)
    z = z
    return theta, r, z

def cart2pol2(x,y,z):
    r = np.sqrt(x**2 + y**2)
    theta = np.atan2(y,x)
    theta = np.unwrap(theta)
    z = z
    return theta, r, z

periodicORaperiodic = 0 #choose 0 if periodic and 1 if aperiodic 
dataPath = '../processedData/'
filePath = '../inputData/'
sSections = 5
Ns = 101 #number o
Nbs = 60
#%%
statorCamber = np.empty([Nbs, sSections, Ns, 3]) #x,y,z
for a in range(Nbs):
    for b in range(sSections):
        statorCamber[a,b,:] = np.loadtxt(dataPath + '/blade{}/blade{}.txt'.format(a+1,b+1), delimiter=',')
                
        
statorCamberCyl = np.empty([Nbs, sSections, Ns, 3]) #theta, r, z
for c in range(Nbs):
    for d in range(sSections):
        for e in range(Ns):
            statorCamberCyl[c,d,e,:] = np.array(cart2pol(statorCamber[c,d,e,0], statorCamber[c,d,e,1], statorCamber[c,d,e,2])).T
        if max(statorCamberCyl[c,d,:,0]) - min(statorCamberCyl[c,d,:,0]) > np.pi:
            statorCamberCyl[c,d,:,:] = np.array(cart2pol2(statorCamber[c,d,:,0], statorCamber[c,d,:,1], statorCamber[c,d,:,2])).T 
  

#%%
chordWisePosStator = np.zeros((Nbs, sSections, Ns, 4))

for f in range(Nbs):
    for g in range(sSections):
        chordWisePosStator[f,g,:,0:3] = statorCamberCyl[f,g,:]
        diffS = np.diff(statorCamberCyl[f,g,:][:,[2,1]], axis=0)        
        seg_lengthS = np.linalg.norm(diffS, axis=1)        
        chordWisePosStator[f,g,:,3] = np.concatenate([[0], np.cumsum(seg_lengthS)])


#%% create list to stack arc lengths
statorChordWisePos = chordWisePosStator.reshape(-1,4)
np.savetxt(filePath + "/xc_stator.txt", statorChordWisePos, delimiter=",")

#%%
# # Nb = 1
# def plot_slicesA(A, B):
#     import matplotlib.pyplot as plt
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     for j in range(Nbs):
#         # j = 0
#         for i in range(sSections):
#             # i = 0
#             ax.plot(A[j,i,:, 0], A[j,i,:, 1], A[j,i,:, 2], 'k', label='tBlade3')
#             # ax.plot(B[j,i,:, 0], B[j,i,:, 1], B[j,i,:, 2], 'r', label='ECL5')
#             # r_A, theta_A, z_A = A[j,i,:,1], A[j,i,:,0], A[j,i,:,2]
#             # x_A, y_A = r_A * np.cos(theta_A), r_A * np.sin(theta_A)
#             # ax.plot(z_A, y_A, x_A, 'k', label='tBlade3')
#             # r_B, theta_B, z_B = B[j,i,:,1], B[j,i,:,0], B[j,i,:,2]
#             # x_B, y_B = r_B * np.cos(theta_B), r_B * np.sin(theta_B)
#             # ax.plot(z_B, y_B, x_B, 'r', label='ECL5')           

#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     plt.axis('equal')
# plot_slicesA(statorCamber, Nbs)    
#%%
# import matplotlib.pyplot as plt

# im = plt.tricontourf(statorChordWisePos[:,2],   # x → axial
#                      statorChordWisePos[:,1],#*statorChordWisePos[:,0],   # y → radial
#                      statorChordWisePos[:,3],   # values → arc length s
#                      cmap='viridis', levels=200)

# # Add colorbar with label in one line
# plt.colorbar(im, label='Arc Length $s$ (m)')

# # Rest of your figure
# plt.axis('equal')
# # plt.xlabel('Axial distance (m)')
# # plt.ylabel('Radial distance (m)')
# plt.xlabel(r'c/$r_{\mathrm{tip}}$')
# plt.ylabel(r'r/$r_{\mathrm{tip}}$')
# plt.title('Rotor camber arc length distribution')

# plt.tight_layout()
# # plt.savefig('./Results/Plots/arcLength.jpg', dpi=1000)
# plt.savefig('/home/adekola/Documents/New_PhD/TestCase4Approach/nonAxiForceFields/plot/arcLength.png', dpi=500, bbox_inches='tight')
# plt.show()