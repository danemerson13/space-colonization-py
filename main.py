import numpy as np
import time

import colony
import kCluster
import plotter

def main():
    # Geometry parameters
    n_clust = 32
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
    init1 = colony.Colony(D, dk, di, in_targets, pts, inlet)
    init1.createTree()
    init1.trimNodes()

    # Build the full tree from the initialization
    col1 = colony.Colony(D, dk, di, kmeans.cluster_centers_, pts, None, init1.nodeList)

    col1.createTree()
    col1.trimNodes()
    col1.createBranches()
    col1.setRadii(initialRadius)
    col1.setResistance(mu)

    # Initialize outlet tree
    init2 = colony.Colony(D, dk, di, out_targets, pts, outlet)
    init2.createTree()
    init2.trimNodes()

    # Build the full tree from the initialization
    col2 = colony.Colony(D, dk, di, kmeans.cluster_centers_, pts, None, init2.nodeList)

    col2.createTree()
    col2.trimNodes()
    col2.createBranches()
    col2.setRadii(initialRadius)
    col2.setResistance(mu)

    # col2.replaceTerminalNodes(col1)

    # col.solveResistanceNetwork(Pin, Pout)

    # plotter.plotFlowRate(col)
    # plotter.plotPressure(col)
    colony.solveResistanceNetwork2x(col1, col2, Pin, Pout)

    plotter.plotFlowRate2x(col1, col2)
    plotter.plotPressure2x(col1, col2)

    print("Hello World")

if __name__ == "__main__":
    main()