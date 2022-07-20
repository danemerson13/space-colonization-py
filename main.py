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
    # n_clusters = np.concatenate([np.arange(5,101,5),np.arange(110,201,10),np.arange(250,501,50)])
    n_clusters = np.arange(110,201,10)
    for n_clust in n_clusters:
        run(n_clust)

def run(n_clust):
    # Geometry parameters
    # n_clust = 100
    initialRadius = 7.5 # mm?

    # Algorithm parameters
    D = 15 # Step distance
    dk = 20 # Kill distance
    di = 25 # Sphere of influence distance

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

    start = time.time()

    # Intialize the first two branches
    init1 = colony.Colony(D, dk, di, in_targets, pts, inlet)
    initTree(init1)
    init2 = colony.Colony(D, dk, di, out_targets, pts, outlet)
    initTree(init2)

    # Build the full tree from the initialization
    col1 = colony.Colony(D, dk, di, kmeans.cluster_centers_, pts, None, init1.nodeList)
    prepTree(col1, initialRadius, mu)
    col2 = colony.Colony(D, dk, di, kmeans.cluster_centers_, pts, None, init2.nodeList)
    prepTree(col2, initialRadius, mu)

    colony.solveResistanceNetwork2x(col1, col2, Pin, Pout)

    end = time.time()

    # plotter.plotFlowRate2x(col1, col2)

    print('Number of Clusters %d, Req = %.3f, V = %.0f, L = %.0f, Time = %.1f secs' %(n_clust, (Pin - Pout)/col1.branchList[0].getFlowRate(), col1.totalVolume(), col1.totalLength(), end - start))

if __name__ == "__main__":
    main()