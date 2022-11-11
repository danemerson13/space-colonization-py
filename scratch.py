import numpy as np
import open3d as o3d

def createVisualizer(objs):
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    for obj in objs:
        vis.add_geometry(obj)
    vis.run()
    vis.destroy_window()
    return vis

def rotate_view(vis):
    ctr = vis.get_view_control()
    ctr.rotate(5.0, 0.0)
    return

def rotateFunc(obj):
    # vis = createVisualizer(obj)
    o3d.visualization.draw_geometries_with_animation_callback([obj], rotate_view)

def main():
    pts = np.load('Point Clouds/liverSample32914.npy')
    pcd = o3d.utility.Vector3dVector(pts)
    pcd = o3d.geometry.PointCloud(pcd)
    rotateFunc(pcd)

if __name__ == "__main__":
    main()