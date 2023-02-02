import numpy as np
import time
import os
import sys
import pickle

from src import colony
from util import plotter, kCluster

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
    nSL = [5,10,25,50,100,250]
    sys.setrecursionlimit(10**6)
    for n in nSL:
        start = time.time()
        path = "data/Point Clouds/10k Clouds/liverSample" + str(0) + ".npy"
        col = createModel(path, nSL = n)
        col.createSegments(lmax = 100)
        col.connectSegments()
        col.fillTree(0.5)
        # Save the conc and times
        tconc = dict(zip(col.tList, col.concentrationList))
        filename = "tconc" + str(n) +".pkl"
        file = open(filename, 'wb')
        pickle.dump(tconc, file)
        end = time.time()
        print("nSL = " + str(n) + ", took " + str(end-start) + " seconds.")

    # for n in nSL:
    #     for i in range(10):
    #         start = time.time()
    #         path = "Point Clouds/liverSample" + str(i) + ".npy"
    #         col = createModel(path, nSL = n)
    #         dir = 'Saved Models' + '/' + 'n' + str(n) + '_' + str(i)
    #         try:
    #             os.mkdir(dir)
    #         except OSError as error:
    #             pass
    #         col.saveModel(dir)
    #         end = time.time()
    #         print("SL = %d, Point Cloud #%d took %.2f secs" %(n, i, end-start))

def main():
    # auto()
    sys.setrecursionlimit(10**6)
    path = "data/Point Clouds/10k Clouds/liverSample" + str(0) + ".npy"
    start = time.time()
    col = createModel(path, nSL = 50)
    end = time.time()
    print("Time to create model: ", end - start)
    # plotter.plotByBranch(col)
    start = time.time()
    col.createSegments(100)
    col.connectSegments()
    end = time.time()
    print("Time to segment: ", end - start)
    start = time.time()
    # col.fillTree(0.5)
    plotter.plotBySegment(col)
    end = time.time()
    print("Time to fill: ", end - start)
    print("All Done!")
    # col.saveModel(os.getcwd())

if __name__ == "__main__":
    main()