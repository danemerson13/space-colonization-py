import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def visualizeClusters(pts,kmeans):
    # check whether 2D or 3D
    n_clusters = np.max(kmeans.labels_)
    if (pts.shape[1] == 2):
        for i in range(n_clusters+1):
            plt.scatter(pts[kmeans.labels_ == i,0], pts[kmeans.labels_ == i,1], marker = 'o', alpha = 0.5, label = 'Cluster ' + str(i))
            plt.scatter(kmeans.cluster_centers_[i,0], kmeans.cluster_centers_[i,1], marker = '*', s = 20, label = "Cluster " + str(i) )
        # plt.legend()
        plt.show()
        # plt.scatter(kmeans.cluster_centers_[i,0], kmeans.cluster_centers_[i,1], marker = 's', size = 2)
    elif (pts.shape[1] == 3):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        for i in range(n_clusters+1):
            ax.scatter(pts[kmeans.labels_ == i,0], pts[kmeans.labels_ == i,1], pts[kmeans.labels_ == i,2], marker = 'o', alpha = 0.5, label = 'Cluster ' + str(i))
            ax.scatter(kmeans.cluster_centers_[i,0], kmeans.cluster_centers_[i,1], kmeans.cluster_centers_[i,2], marker = '*', s = 20)
        ax.set_box_aspect((np.ptp(pts[:,0]), np.ptp(pts[:,1]), np.ptp(pts[:,2])))
        ax.view_init(-170,-75)
        # plt.legend()
        plt.show()
    else:
        raise TypeError("Incorrect Dimension")

def myKMeans(pts, n_clust):
    cluster_percent = 0.05
    # n_clust = 100
    n = pts.shape[0]
    kmeans = KMeans(n_clusters = n_clust, max_iter = 1000, random_state = 0, n_init = 10).fit(pts)
    # kmeans = KMeans(n_clusters = int(pts.shape[0]//(n*cluster_percent)), max_iter = 1000, random_state = 0).fit(pts)
    return kmeans

def main():
    pts = np.load('Point Clouds/liverSample3286.npy')
    n_clust = 8
    kmeans = myKMeans(pts, n_clust)
    visualizeClusters(pts,kmeans)

if __name__ == "__main__":
    main()