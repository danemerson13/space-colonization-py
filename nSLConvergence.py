import numpy as np
import time
import os
import sys
import pickle

from src import colony
from util import plotter
from util.KMeans import kCluster

'''
In this script we will perform a convergence study on the number of SLs
nSLs = [5, 10, 25, 50, 100, 250, 500, 1000, 2000]
We fix a time step based on the tConvergence study, dt = 0.01s

Then we need to generate a set of models from varying nSLs
The point cloud should be fixed, and we will use a dense point cloud

We will want to record the following:
Time to create model (tSpaceCol)
Time to solve matrix once (tMatrixSolve)
Time to solve pressures/flows (tRSLSolve)
NOTE: the time for these two is affect heavily by the fact the matrix solve is implemented as a dense matrix rather than a sparse one

All of the below can be solved after model generation (and resimulated for different time steps, dt)
Time to simulate filling (tFillSolve)
Time to 50% total concentration
Time to 98% total concentration
Time for last SL > 50% concentration
Time for last SL > 98% concentration
Total concentration at t = 200s
'''

def createModel(pointPath, nSL):
    # Geometry Parameters
    n_clust = nSL
    initialRadius = 7.5 # mm

    # Algorithm Parameters
    D = 5 # Step distance
    dk = 10 # Kill distance
    di = 15 # Sphere of influence distance

    # Fluid parameters
    mu = 3.5 # centipoise
    rho = 1 # g/mL

    # Model Parameters
    Pin = 5.8 # mmHg
    Pout = 0 # mmHg
    Qin = 767 # mL/min

    # Import the uniformly sampled point cloud on the liver volume
    # This will be used as the attractors and also to generate the point clouds
    pts = np.load(pointPath)

    # Define inlet and outlet points to serve as root nodes for the two trees
    inlet = np.array([-8.15, 18.59, 159.07])
    outlet = np.array([-7.82, 4.37, 3.47])

    # Define additional target points for realistic initializion
    inlet_targets = np.array([[-50,0,110],[40,20,90]])
    outlet_targets = np.array([[-60,-10,75],[-10,30,75],[50,15,75]])

    # Call KMeans to get the Super Lobule locations
    kmeans = kCluster.myKMeans(pts, n_clust)

    # Initialize trees
    init1 = colony.Colony(D,dk,di)
    init1.initTree(inlet_targets, pts, inlet)
    init2 = colony.Colony(D,dk,di)
    init2.initTree(outlet_targets, pts, outlet)

    # Create full trees from initializations
    col1 = colony.Colony(D,dk,di)
    col1.createTree(kmeans.cluster_centers_, pts, init1.nodeList, initialRadius, mu)
    col2 = colony.Colony(D,dk,di)
    col2.createTree(kmeans.cluster_centers_, pts, init2.nodeList, initialRadius, mu)

    # Merge trees
    col1.mergeTrees(col2)

    # Solve
    col1.solveRSL(Qin, Pin, Pout, n_iter=1000, tol=1e-6, verbosity=0, a=0, b=1e-4)
    col1.setSLVolume(375000) # Setting SL volume with bloodVol = 375 mL

    return col1

def auto():
    nSL = [5, 10, 25, 50, 100, 250, 500, 1000, 2000]
    sys.setrecursionlimit(10**9)
    for n in nSL:
        print("Creating %dSL Model" %(n))
        # Use the fixed rich point cloud for consistency
        pointpath = 'data/Point Clouds/100k Clouds/liverSample' + str(0) + '.npy'
        # Create the model, simulate filling, save
        start = time.time()
        col = createModel(pointpath, nSL = n)
        end = time.time()
        print("Model created in ", end - start, " seconds.")
        start = time.time()
        col.fillTree(dt = 0.01, tStop = 200, gif = False)
        end = time.time()
        print("Tree filled in ", end - start, " seconds.")
        start = time.time()
        folder = os.getcwd() + '/results/nSLConvergence/' + str(n) + 'SL'
        # Check that directory exists, if not make it
        if not os.path.exists(folder):
            os.mkdir(folder)
        col.saveModel(path = folder)
        end = time.time()
        print("Model saved in ", end - start, " seconds.")

if __name__ == "__main__":
    auto()