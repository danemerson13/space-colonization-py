import numpy as np
import time

import colony
import kCluster
import plotter

def initTree(col):
    col.createTree(initflag = True)
    col.trimNodes()

def prepTree(col, initialRadius, mu):
    col.createTree(initflag = False)
    col.trimNodes()
    col.createBranches()
    col.setRadii(initialRadius)
    col.setResistance(mu)

def main():
    # Geometry parameters
    n_clust = 10
    initialRadius = 7.5 # mm?

    # Algorithm parameters
    D = 10 # Step distance
    dk = 15 # Kill distance
    di = 20 # Sphere of influence distance

    # Blood parameters
    mu = 3.5 # centipoise
    Pin = 10 # mmHg
    Pout = 0 # mmHg

    # Import the uniformly sampled point cloud on the liver volume
    # This will be used as the attractors and also to generate the point clouds
    pts = np.load('Point Clouds/liverSample3286.npy')

    # Define Inlet and Outlet points to serve as the root nodes
    inlet = np.array([-5.52, 3.08, 2.45])
    outlet = np.array([-13.49, 14.48, 114.21])

    # Call KMeans to get the Super Lobule locations
    kmeans = kCluster.myKMeans(pts, n_clust)

    # Inlet targets 
    in_targets = np.array([[-50., 0., 50.], [35., 0., 45.]])
    # Outlet targets
    out_targets = np.array([[25., 20., 60.], [-35., 10., 80.]])

    # Intialize the first two branches
    init1 = colony.Colony(D, dk, di, "Inlet", in_targets, pts, inlet)
    initTree(init1)

    init2 = colony.Colony(D, dk, di, "Inlet", out_targets, pts, outlet)
    initTree(init2)

    # plotter.plotByNode(init1)
    # plotter.plotByNode(init2)

    # Build the full tree from the initialization
    col1 = colony.Colony(D, dk, di, "Inlet", kmeans.cluster_centers_, pts, None, init1.nodeList)
    prepTree(col1, initialRadius, mu)

    col2 = colony.Colony(D, dk, di, "Inlet", kmeans.cluster_centers_, pts, None, init2.nodeList)
    prepTree(col2, initialRadius, mu)

    out = colony.reverseColony(col2)

    plotter.plotByBranch(col2)
    plotter.plotByBranch(out)

    print("Hello World")

if __name__ == "__main__":
    main()