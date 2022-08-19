import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from stl import mesh as np_mesh
from mpl_toolkits import mplot3d

def plotByBranch(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'k-')

        if branch.getType() == "Inlet":
            if branch.isSL():
                ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'bo')
            elif branch.isRoot():
                ax.plot3D(coords[0][0],coords[0][1],coords[0][2],'go')
        elif branch.getType() == "Outlet":
            if branch.isSL():
                ax.plot3D(coords[0][0],coords[0][1],coords[0][2],'bo')
            elif branch.isTerminal():
                ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'ro')

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    plt.show()

def plotPressure(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min, max = getPressureLims(col)
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'k-')

    for node in col.nodeList:
        if node.getType() != "Interior":
            nodeLoc = node.getLocation()
            ax.plot3D(nodeLoc[0], nodeLoc[1], nodeLoc[2], 'o', color = mymap.to_rgba(node.getPressure())[:-1])

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Pressure', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.show()

def plotFlowRate(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min, max = getFlowRateLims(col)
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        color = mymap.to_rgba(branch.getFlowRate())[:-1]
        rad = branch.getRadius()
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Flow Rate', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.show()

##### HELPER FUNCTIONS #####

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

def getPressureLims(col):
    max = -np.inf
    min = np.inf
    for node in col.nodeList:
        if node.getType() != "Interior":
            PNode = node.getPressure()
            if PNode > max:
                max = PNode
            if PNode < min:
                min = PNode
    return min, max

def getFlowRateLims(col):
    max = -np.inf
    min = np.inf
    for branch in col.branchList:
        QBranch = branch.getFlowRate()
        if QBranch > max:
            max = QBranch
        if QBranch < min:
            min = QBranch
    return min, max