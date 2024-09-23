import numpy as np

class CubicBezierCurve:
    def __init__(self, p0, p1, p2, p3):
        self.P0 = p0
        self.P1 = p1
        self.P2 = p2
        self.P3 = p3

    def eval(self, parameterization, tangent = False):
        # Parameterization must be an array of values ∈ [0,1]
        trajectory = np.zeros((len(parameterization), len(self.P0)))
        if tangent:
            tangents = np.zeros((len(parameterization), len(self.P0)))
        for i, t in enumerate(parameterization):
            # Use formula for cubic Bezier curve
            trajectory[i,:] = (1 - t)**3*self.P0 + 3*(1 - t)**2*t*self.P1 + 3*(1-t)*t**2*self.P2 + t**3*self.P3
            if tangent:
                # tangents[i,:] = -3*(1-t)**2*self.P0 + 3**(1-t)**2*self.P1 - 6*t*(1-t)*self.P1 -3*t**2*self.P2 + 6*t*(1-t)*self.P2 + 3*t**2*self.P3
                tangents[i,:] = -3*(t-1)**2*self.P0 + 3*(t-1)**2*self.P1 + 6*(t-1)*t*self.P1 - 3*t**2*self.P2 - 6*(t-1)*t*self.P2 + 3*t**2*self.P3

        if tangent:
            return trajectory, tangents
        else:
            return trajectory
    
class G1CubicSpline:
    def __new__(self, points):
        # For a set of points, define a piecewise series of bezier curves interpolating with G1 continuity
        # Points is a list of points to be interpolated
        num_pts = len(points)
        bezierList = list()
        # If less than 2 points are passed raise an error
        if num_pts < 2:
            raise ValueError("Number of points must be at least 2")
        # If there are just two points make a straight line
        elif num_pts == 2:
            p0 = points[0]
            p1 = points[0] * 2/3 + points[1] * 1/3
            p2 = points[0] * 1/3 + points[1] * 2/3
            p3 = points[1]
            bezierList.append(CubicBezierCurve(p0,p1,p2,p3))
            return bezierList
        # Perform G1 piecewise interpolation for any other number of points
        else:
            # Extract the individual x and y coordinates
            xcoord = [point[0] for point in points]
            ycoord = [point[1] for point in points]
            zcoord = [point[2] for point in points]
            # Compute the Bezier coefficients
            cx = fitBezier(xcoord, num_pts-1)
            cy = fitBezier(ycoord, num_pts-1)
            cz = fitBezier(zcoord, num_pts-1)
            # Create a list of the Bezier coeffs
            coeffs = list()
            for i in range(num_pts):
                coeffs.append([cx[i], cy[i], cz[i]])
            # Finally compute the control points of each individual Bezier curve
            for i in range(num_pts-1):
                p0 = points[i+0]
                p1 = points[i+0] + coeffs[i+0]
                p2 = points[i+1] - coeffs[i+1]
                p3 = points[i+1]
                bezierList.append(CubicBezierCurve(p0,p1,p2,p3))
            return bezierList
    
def fitBezier(a, n):
    # Given the points in a, determine the bezier coefficients with G1 continuity
    d1 = np.zeros(n+1) # tangent vector at each control point
    d2 = np.zeros(n+1) # second derivative at eaach control point
    coeff = np.zeros(n+1) # coefficients of each control point
    
    # Compute the first derivatives
    d1[0] = 2 * (a[1] - a[0]) / 3
    for i in range(1,n):
        d1[i] = a[i+1] - a[i-1]
    d1[n] = 2 * (a[n] - a[n-1]) / 3

    # Compute the second derivatives
    d2[n] = 1
    d2[n-1] = 1
    for i in range(n-2, -1, -1):
        d2[i] = 4 * d2[i+1] - d2[i+2]
    for i in range(1, n+1, 2):
        d2[i] *= -1

    # Compute the coefficients
    coeff[0] = 0
    for i in range(n, -1, -1):
        coeff[0] += d1[i] * d2[i]
    coeff[0] /= (d2[0] + d2[1])
    coeff[1] = d1[0] - coeff[0]
    for i in range(1, n):
        coeff[i+1] = d1[i] - 4 * coeff[i] - coeff[i-1]

    return coeff