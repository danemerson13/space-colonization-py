import numpy as np
import open3d as o3d
import matplotlib.pyplot as plt
import matplotlib
from stl import mesh as np_mesh
from mpl_toolkits import mplot3d

def vertexSelect():
    mesh = o3d.io.read_triangle_mesh('Models/liver159.stl')
    mesh.compute_vertex_normals()
    o3d.visualization.draw_geometries_with_vertex_selection([mesh])

def getLiverSTL():
    # Load the STL files and add the vectors to the plot
    mesh = np_mesh.Mesh.from_file('Models/liver159.stl')
    liver = mplot3d.art3d.Poly3DCollection(mesh.vectors)
    liver.set_alpha(0.05)
    liver.set_facecolor([0.5, 0.5, 0.5, 0.25])
    liver.set_edgecolor([0., 0., 0., 0.25])

    return liver, mesh

def bound_mesh(ax, mesh):
    pts = mesh.points.reshape(len(mesh.points)*3,3)

    x_range = abs(max(pts[:,0]) - min(pts[:,0]))
    x_middle = np.mean([max(pts[:,0]), min(pts[:,0])])
    y_range = abs(max(pts[:,1]) - min(pts[:,1]))
    y_middle = np.mean([max(pts[:,1]), min(pts[:,1])])
    z_range = abs(max(pts[:,2]) - min(pts[:,2]))
    z_middle = np.mean([max(pts[:,2]), min(pts[:,2])])

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def visualizeInit():
    inlet = np.array([-8.15, 18.59, 159.07])
    inlet_init = list([np.array([-50,0,110]),np.array([40,20,90])])
    outlet = np.array([-7.82, 4.37, 3.47])
    outlet_init = list([np.array([-60,-10,75]),np.array([-10,30,75]),np.array([50,15,75])])

    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    ax.plot3D(inlet[0], inlet[1], inlet[2], 'go')
    ax.plot3D(outlet[0], outlet[1], outlet[2], 'ro')
    for init in outlet_init:
        ax.plot3D(init[0], init[1], init[2], 'bo')

    bound_mesh(ax, mesh)
    ax.view_init(elev = -165., azim = -120.)
    # ax.view_init(elev=25., azim=-55.)
    plt.show()

def main():
    # vertexSelect()
    visualizeInit()

if __name__ == "__main__":
    main()