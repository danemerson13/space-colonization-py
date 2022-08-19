import numpy as np
import time

import colony
import kCluster
import plotter

def main():
    # Geometry parameters
    n_clust = 4
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
    # pts =np.load('Point Clouds/liverSample32914.npy')

    # Define Inlet and Outlet points to serve as the root nodes
    inlet = np.array([-8.15, 18.59, 159.07])
    outlet = np.array([-7.82, 4.37, 3.47])

    inlet_targets = np.array([[-50,0,110],[40,20,90]])
    outlet_targets = np.array([[-60,-10,75],[-10,30,75],[50,15,75]])

    # Call KMeans to get the Super Lobule locations
    kmeans = kCluster.myKMeans(pts, n_clust)

    # Initialize and construct the inlet and outlet trees
    init1 = colony.Colony(D,dk,di)
    init1.initTree(inlet_targets, pts, inlet)
    col1 = colony.Colony(D,dk,di)
    col1.createTree(kmeans.cluster_centers_, pts, init1.nodeList, initialRadius, mu)

    init2 = colony.Colony(D,dk,di)
    init2.initTree(outlet_targets, pts, outlet)
    col2 = colony.Colony(D,dk,di)
    col2.createTree(kmeans.cluster_centers_, pts, init2.nodeList, initialRadius, mu)

    # Merge the two trees
    col1.mergeTrees(col2)

    col1.solveResistanceNetwork(Pin, Pout)

    plotter.plotPressure(col1)
    plotter.plotFlowRate(col1)

    print("All Done")

if __name__ == "__main__":
    main()