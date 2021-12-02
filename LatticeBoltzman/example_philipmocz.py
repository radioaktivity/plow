import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import glob
import matplotlib.cm as cm
"""
Create Your Own Lattice Boltzmann Simulation (With Python)
Philip Mocz (2020) Princeton Univeristy, @PMocz

Simulate flow past cylinder
for an isothermal fluid

"""

def generate_video(name_video):
    os.chdir("LatticeBoltzman/video")
    r'ffmpeg -framerate 10 -pattern_type glob -i "file*.png" output.mp4'
    subprocess.call([
        'ffmpeg', '-framerate', '25', '-pattern_type', 'glob', '-i', 'file*.png', name_video
    ])
    for file_name in glob.glob("*.png"):
        os.remove(file_name)
        

def main(iterations=5000, write_intervall=10):
    """ Finite Volume simulation """

    # Simulation parameters
    Nx                     = 600   # resolution x-dir
    Ny                     = 600   # resolution y-dir
    rho0                   = 100    # average density
    tau                    = 0.6    # collision timescale
    Nt                     = iterations   # number of timesteps
    plotRealTime = True # switch on for plotting as the simulation goes along

    # Lattice speeds / weights
    NL = 9
    idxs = np.arange(NL)
    cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
    cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
    weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) # sums to 1

    # Initial Conditions
    F = np.ones((Ny,Nx,NL)) #* rho0 / NL
    np.random.seed(42)
    F += 0.01*np.random.randn(Ny,Nx,NL)
    X, Y = np.meshgrid(range(Nx), range(Ny))
    F[:,:,3] += 2 * (1+0.2*np.cos(2*np.pi*X/Nx*4))
    rho = np.sum(F,2)
    for i in idxs:
        F[:,:,i] *= rho0 / rho

    # Cylinder boundary
    X, Y = np.meshgrid(range(Nx), range(Ny))
    cylinder = np.zeros(X.shape)
    for i in range(cylinder.shape[0]):
        for j in range(cylinder.shape[1]):
            if ((i>Nx/2-50) and (i<Nx/2+50)) and ((j>Ny/2-50) and (j<Ny/2+50)):
                cylinder[i,j] = 1

    cylinder = cylinder == 1


    # Prep figure
    fig = plt.figure(figsize=(4,2), dpi=80)

    # Simulation Main Loop
    for it in range(Nt):
        # Drift
        for i, cx, cy in zip(idxs, cxs, cys):
            F[:,:,i] = np.roll(F[:,:,i], cx, axis=1)
            F[:,:,i] = np.roll(F[:,:,i], cy, axis=0)
        
        
        # Set reflective boundaries
        bndryF = F[cylinder,:]
        bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]

        
        # Calculate fluid variables
        rho = np.sum(F,2)
        ux  = np.sum(F*cxs,2) / rho
        uy  = np.sum(F*cys,2) / rho
        
        
        # Apply Collision
        Feq = np.zeros(F.shape)
        for i, cx, cy, w in zip(idxs, cxs, cys, weights):
            Feq[:,:,i] = rho * w * ( 1 + 3*(cx*ux+cy*uy)  + 9*(cx*ux+cy*uy)**2/2 - 3*(ux**2+uy**2)/2 )
        
        F += -(1.0/tau) * (F - Feq)
        
        # Apply boundary 
        F[cylinder,:] = bndryF
        
        
        # plot in real time - color 1/2 particles blue, other half red
        if (plotRealTime and (it % write_intervall) == 0) or (it == Nt-1):
            if True:
                plt.cla()
                ux[cylinder] = 0
                uy[cylinder] = 0
                cmap = plt.cm.bwr
                cmap.set_bad('black')
                umean = np.sqrt(np.multiply(ux,ux) + np.multiply(uy,uy))
                vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) - (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
                vorticity[cylinder] = np.nan
                plt.imshow(rho)#, cmap='bwr')
                plt.clim(-.1, .1)
                ax = plt.gca()
                ax.invert_yaxis()
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)	
                ax.set_aspect('equal')
                plt.savefig("LatticeBoltzman/video/file%02d.png" % it, dpi=240)

            if False:
                plt.cla()
                ux[cylinder] = 0
                uy[cylinder] = 0
                vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) - (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
                vorticity[cylinder] = np.nan
                cmap = plt.cm.bwr
                cmap.set_bad('black')
                plt.imshow(vorticity, cmap='bwr')
                plt.clim(-.1, .1)
                ax = plt.gca()
                ax.invert_yaxis()
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)	
                ax.set_aspect('equal')	
                plt.pause(0.001)
        
        #plt.show()
    return 0


if __name__== "__main__":
    main(iterations=9999, write_intervall=50)
    generate_video("new_animation003.mp4")
