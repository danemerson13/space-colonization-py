from matplotlib.pyplot import plot
import numpy as np
import time

import colony
import kCluster
import plotter
import scipy

# def main():
# # # # #     # n_clusters = list([5,10,15,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300])
# # # # #     n_clusters = list([5,10,15,25,40,60,100,150,200,300])
#     n_clusters = list([5,10,25,50,100,250])
#     for n_clust in n_clusters:
#         run(n_clust)

# def run(n_clust):
def main():
    # Geometry parameters
    n_clust = 25
    initialRadius = 7.5 # mm?

    # Algorithm parameters
    D = 5 # Step distance
    dk = 10 # Kill distance
    di = 15 # Sphere of influence distance

    # Blood parameters
    mu = 3.5 # centipoise
    rho = 1 # g/mL
    Pin = 5.8 # mmHg
    Pout = 0 # mmHg
    Qin = 767 # mL/min

    # Import the uniformly sampled point cloud on the liver volume
    # This will be used as the attractors and also to generate the point clouds
    pts = np.load('Point Clouds/liverSample4.npy')
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

    col1.solveRSL(Qin, Pin, Pout, n_iter = 1000, tol = 1e-6, verbal = 0)

    # col1.solveResistanceNetwork(Pin, Pout, flag = True)

    plotter.plotFlowRate(col1)

    # col1.overfillNetwork(0.1)

    # print("SL: ", n_clust, ", Vol:", col1.totalVolume())
    # col1.setSLVolume(375000) # volume in mm3

    # mergedCounter = col1.setGeneration(0.9, 120*np.pi/180)
    # col1.generationStatistics(0.7, 120)
    # plotter.plotGenerations(col1)

    # col1.generationBoxPlots(0.7,120)

if __name__ == "__main__":
    main()