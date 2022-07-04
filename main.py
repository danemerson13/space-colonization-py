import numpy as np
import time

import colony
import kCluster
import plotter

def main():
    # Geometry parameters
    n_clust = 100
    initialRadius = 7.5 # mm?

    # Algorithm parameters
    D = 5 # Step distance
    dk = 10 # Kill distance
    di = 15 # Sphere of influence distance

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
    init = colony.Colony(D, dk, di, in_targets, pts, inlet)
    init.createTree()
    init.trimNodes()

    # Build the full tree from the initialization
    col = colony.Colony(D, dk, di, kmeans.cluster_centers_, pts, None, init.nodeList)

    col.createTree()
    col.trimNodes()
    col.createBranches()
    col.setRadii(initialRadius)
    col.setResistance(mu)

    col.solveResistanceNetwork(Pin, Pout)

    # plotter.plotFlowRate(col)
    # plotter.plotPressure(col)

    # col.assignGenerations()
    col.assignGenerationsThreshold(0.65)

    col.computeGeometryStatistics()

    plotter.plotGeneration(col)

    print("Hello World")

if __name__ == "__main__":
    main()