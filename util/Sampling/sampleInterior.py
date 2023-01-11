import open3d as o3d
import numpy as np
import myRandom
import mollerTrumbore
import time

def loadMesh(path):
    mesh = o3d.io.read_triangle_mesh(path)
    return mesh

def visualizeMesh(obj):
    if type(obj) == o3d.cpu.pybind.geometry.TriangleMesh:
        obj.compute_vertex_normals()
    o3d.visualization.draw_geometries([obj])

def samplePts(mesh,n):
    # Get bounding box
    min = mesh.get_min_bound()
    max = mesh.get_max_bound()
    rg = max - min
    pts = myRandom.random(3,n) * rg + min
    return pts

def randomRay():
    ray = np.random.rand(3)-0.5
    return ray/np.linalg.norm(ray)

def getTriangleCoords(mesh):
    T = np.asarray(mesh.triangles)
    V = np.asarray(mesh.vertices)
    return V[T]

def countIntersections(mesh,point):
    triangles = getTriangleCoords(mesh)
    ray = randomRay()
    count = 0
    for i in range(len(triangles)):
        count += mollerTrumbore.rayTriangleIntersect(point, ray, triangles[i])
    return count

def checkPts(mesh,pts):
    mask = np.zeros(len(pts), dtype = bool)
    for i in range(len(pts)):
        # if not i % 1000:
        #     print(i , "points checked")
        count = countIntersections(mesh,pts[i])
        # if odd num of intersections, pt is inside mesh
        if not(count % 2 == 0):
            mask[i] = True
    return pts[mask]

def main(n):
    path = '../../Models/liver159.stl'
    mesh = loadMesh(path)
    # visualizeMesh(mesh)
    # n = 10000
    print('Sampling ' + str(n) + ' points...')
    pts = samplePts(mesh,n)
    start = time.time()
    pts = checkPts(mesh, pts)
    end = time.time()
    print("It took ", end - start, " seconds to check the points")
    # A = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(pts))
    # visualizeMesh(A)
    name = 'liverSample' + str(len(pts))
    np.save(name, pts)
    print('All Done')

def automatic(samples, npts):
    path = '../../Models/liver159.stl'
    mesh = loadMesh(path)
    for i in range(samples):
        start = time.time()
        pts = samplePts(mesh,npts)
        pts = checkPts(mesh, pts)
        name = 'liverSample' + str(i)
        np.save(name, pts)
        end = time.time()
        print("It took ", end-start, " seconds to sample and check the points for point cloud #", i, ". ", len(pts), " points are in the cloud.")

if __name__ == "__main__":
    main()
    # automatic(10,10000)