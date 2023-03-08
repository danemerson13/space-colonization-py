import numpy as np
import time
import os
import sys
import pickle

from src import colony
from util import plotter
from util.KMeans import kCluster

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
    nSL = [5,10,25,50,100,250,500,1000,2000]
    sys.setrecursionlimit(10**9)
    for n in nSL:
        # Check that directory exists, if not make it
        if not os.path.exists(os.getcwd() + '/results/' + str(n) + 'SL'):
            os.mkdir(os.getcwd() + '/results/' + str(n) + 'SL')
        for i in range(10):
            print("%dSL%d" %(n, i))
            # Use richer point cloud for models with more SLs
            # if n <= 100:
            #     path = 'data/Point Clouds/10k Clouds/liverSample' + str(i) + '.npy'
            # else:
            path = 'data/Point Clouds/100k Clouds/liverSample' + str(i) + '.npy'
            # Create the model, simulate filling, save
            start = time.time()
            col = createModel(path, nSL = n)
            end = time.time()
            print("Model created in ", end - start, " seconds.")
            start = time.time()
            col.fillTree(dt = 0.1)
            end = time.time()
            print("Tree filled in ", end - start, " seconds.")
            start = time.time()
            folder = os.getcwd() + '/results/' + str(n) + 'SL' + '/sample' + str(i)
            # Check that directory exists, if not make it
            if not os.path.exists(folder):
                os.mkdir(folder)
            col.saveModel(path = folder)
            end = time.time()
            print("Model saved in ", end - start, " seconds.")

def main():
    i = 0 
    n = 250

    # Check that directory exists, if not make it
    # if not os.path.exists(os.getcwd() + '/trial/' + str(n) + 'SL'):
    #     os.mkdir(os.getcwd() + '/trial/' + str(n) + 'SL')

    sys.setrecursionlimit(10**9)

    print("%dSL%d" %(n, i))
    # Use richer point cloud for models with more SLs
    if n <= 100:
        path = 'data/Point Clouds/10k Clouds/liverSample' + str(i) + '.npy'
    else:
        path = 'data/Point Clouds/100k Clouds/liverSample' + str(i) + '.npy'
    # Create the model, segment, simulate filling, save model
    start = time.time()
    col = createModel(path, nSL = n)
    end = time.time()
    print("Model created in ", end - start, " seconds.")
    start = time.time()
    col.fillTree(dt = 0.2, gif = True)
    end = time.time()
    print("Filling algorithm took ", end- start, "seconds.")
    # Save the tree
    start = time.time()
    folder = os.getcwd() + '/misc/trial/' + str(n) + 'SL' + '/sample' + str(i)
    # Check that directory exists, if not make it
    if not os.path.exists(folder):
        os.mkdir(folder)
    col.saveModel(path = folder)
    end = time.time()
    print("Model saved in ", end - start, " seconds.")

if __name__ == "__main__":
    # main()
    auto()