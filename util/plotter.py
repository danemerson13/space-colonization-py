import open3d as o3d
import numpy as np
import copy
import matplotlib
import matplotlib.pyplot as plt
from stl import mesh as np_mesh
from mpl_toolkits import mplot3d
import sys
sys.path.append('../../Space Colonization/')
from src import superlobule, branch
import os

def plotConcentration(col, time = None, path = None):
    # Create the figure
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for obj in col.fillList:
        # Interpolate color
        conc = obj.getConcentration()
        color = conc*np.array([1,0,0]) + (1-conc)*np.array([0,0,1])
        if isinstance(obj, branch.Branch):
            coords = getBranchCoords(obj)
            rad = obj.getRadius()
            color = conc*np.array([1,0,0]) + (1-conc)*np.array([0,0,1])
            for j in range(len(coords) - 1):
                ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)
        else: # isinstance(obj, superlobule.SuperLobule)
            loc = obj.getLocation()
            ax.plot3D(loc[0],loc[1],loc[2],'o', ms = 2, color = color)

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.text(25,0,150,'Concentration: %.3f' %(conc))
    ax.text(25,0,175,"Time: %.2f s" %(time))
    ax.set_axis_off()
    if path:
        plt.savefig(path, bbox_inches = 'tight', dpi = 100)
        plt.close()
    else:
        plt.show()

def plotFlowRate(col, path = None):
    # Create color mapping, max note to convert flowRate from mm^3/s to mmHg/min by multiplying by 3/50
    fmin, fmax = getFlowRateLims(col)
    cmap = matplotlib.cm.viridis
    normalize = matplotlib.colors.Normalize(vmin = fmin*3/50, vmax = fmax*3/50)
    mymap = matplotlib.cm.ScalarMappable(norm = normalize, cmap = cmap)
    
    # Create the figure
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        color = mymap.to_rgba(branch.getFlowRate() * 3/50)[:-1]
        rad = branch.getRadius()
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Flow Rate (mL/min)', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    if path:
        plt.savefig(path, dpi = 200)
        plt.close()
    else:
        plt.show()

def plotBranch3D(col):
    meshList = list()
    for branch in col.branchList:
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            if branch.getType() == "Inlet":
                meshList.append(createCylinder(coords[j],coords[j+1],branch.getRadius(), np.array([1,0,0])))
            if branch.getType() == "Outlet":
                meshList.append(createCylinder(coords[j],coords[j+1],branch.getRadius(), np.array([0,0,1])))

    o3d.visualization.draw_geometries(meshList)

##### Helper Functions #####

def createConeAdv(baseCenter, topCenter, baseRadius, topRadius, baseAxis, topAxis, stackCount, sliceCount):
    # Functionally the same as createCone() but also accounts for the normal direction 
    # of each cone face and interpolates intermediate stack normals

    # Check slice and stack count is 2 or more
    if sliceCount < 2:
        raise ValueError("Slice count must be 2 or higher")
    if stackCount < 2:
        raise ValueError("Stack count must be 2 or higher")

    V = list()
    F = list()

    axis = topCenter - baseCenter
    coneLength = np.linalg.norm(axis)
    axis = axis/coneLength
    
    dTheta = 2 * np.pi / sliceCount
    # Normalize base and top axes
    baseAxis = baseAxis/np.linalg.norm(baseAxis)
    topAxis = topAxis/np.linalg.norm(topAxis)
    # Make sure random ray is orthogonal to the two axes
    randomRay = np.cross(baseAxis,topAxis)

    # Make the vertices
    for i in range(stackCount):
        # Compute center of the stack
        center = baseCenter + i/stackCount * coneLength * axis
        # Interpolate the radius based on stack
        stackRadius = baseRadius + (topRadius - baseRadius) * i/(stackCount-1)
        print(stackRadius)
        # Determine stack axis by interpolating between base and stack axes
        stackAxis = (1 - i/(stackCount-1)) * baseAxis + i/(stackCount-1) * topAxis
        stackAxis = stackAxis/np.linalg.norm(stackAxis)
        print(stackAxis)
        # Rotate around stack generating vertex for each slice
        for j in range(sliceCount):
            # TODO inconsistent rotation of face because of np.cross(np.cross())
            V.append(center + stackRadius * (np.cos(dTheta*j)*np.cross(stackAxis,randomRay) + np.sin(dTheta*j)*np.cross(np.cross(stackAxis,randomRay),axis)))
            # Assign the indices
            if i < stackCount-1 and j < sliceCount-1:
                # F.append([sliceCount*i+j, sliceCount*i+j+1, sliceCount*(i+1)+j+1])
                # F.append([sliceCount*i+j, sliceCount*(i+1)+j+1, sliceCount*(i+1)+j])
                F.append([sliceCount*i+j, sliceCount*(i+1)+j+1, sliceCount*i+j+1])
                F.append([sliceCount*i+j, sliceCount*(i+1)+j, sliceCount*(i+1)+j+1])
        if i < stackCount - 1:
            # F.append([sliceCount*(i+1)-1, sliceCount*i, sliceCount*(i+1)])
            # F.append([sliceCount*(i+1)-1, sliceCount*(i+1), sliceCount*(i+2)-1])
            F.append([sliceCount*(i+1)-1, sliceCount*(i+1), sliceCount*i])
            F.append([sliceCount*(i+1)-1, sliceCount*(i+2)-1, sliceCount*(i+1)])

    # Make the mesh
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(V)
    mesh.triangles = o3d.utility.Vector3iVector(F)
    mesh.compute_triangle_normals()
    mesh.compute_vertex_normals()

    return mesh

def createSphere(loc, radius, color):
    # Make template sphere
    sph = o3d.geometry.TriangleMesh.create_sphere(radius)
    # Translarte to correct position
    sph = copy.deepcopy(sph).translate(loc)
    # Make it look nice
    sph.paint_uniform_color(color)
    sph.compute_triangle_normals()
    sph.compute_vertex_normals()
    return sph

def createCylinder(a, b, radius, color):
    # Make template cylinder aligned with z-axis
    height = np.linalg.norm(b-a)
    cyl = o3d.geometry.TriangleMesh.create_cylinder(radius, height)
    # Need to rotate and translate the cylinder to its correct location
    # Compute axis and angle between vectors
    vec1 = np.array([0,0,1])
    vec2 = (b - a) / np.linalg.norm(b-a)
    axis = np.cross(vec1, vec2) / np.linalg.norm(np.cross(vec1, vec2))
    angle = np.arccos((np.dot(vec1, vec2)))
    # Apply rotation and translation
    R = cyl.get_rotation_matrix_from_axis_angle(axis * angle)
    cyl = copy.deepcopy(cyl).rotate(R, center = (0,0,0))
    cyl = copy.deepcopy(cyl).translate((a+b)/2)
    # Make it look nice
    cyl.paint_uniform_color(color)
    cyl.compute_triangle_normals()
    cyl.compute_vertex_normals()
    return cyl

def getBranchCoords(branch):
    coordList = list()
    # Check the type of branch. 
    # For Inlet branches go from distal to proximal using getParents()[0]
    # For Outlet branches go from proximal to distal using getChildren()[0]
    if branch.getType() == "Inlet":
        curNode = branch.getDistal()
        coordList.append(curNode.getLocation())
        while curNode != branch.getProximal():
            curNode = curNode.getParents()[0]
            coordList.append(curNode.getLocation())
        # For the inlet branches reverse the list so we go from root -> terminal
        coordList.reverse()
    elif branch.getType() == "Outlet":
        curNode = branch.getProximal()
        coordList.append(curNode.getLocation())
        while curNode != branch.getDistal():
            curNode = curNode.getChildren()[0]
            coordList.append(curNode.getLocation())
    else:
        raise TypeError("Branch has no type")
    return coordList

def getLiverSTL():
    # Load the STL files and add the vectors to the plot
    mesh = np_mesh.Mesh.from_file(os.path.abspath(os.path.join(os.getcwd(), '..'))+ '/data/Models/liver159.stl')
    liver = mplot3d.art3d.Poly3DCollection(mesh.vectors)
    liver.set_alpha(0.05)
    liver.set_facecolor([0.5, 0.5, 0.5, 0.20])
    liver.set_edgecolor([0., 0., 0., 0.20])

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

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
    ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def getFlowRateLims(col):
    maxFlow = -np.inf
    minFlow = np.inf
    for branch in col.branchList:
        flow = branch.getFlowRate()
        if flow > maxFlow:
            maxFlow = flow
        if flow < minFlow:
            minFlow = flow
    return minFlow, maxFlow
