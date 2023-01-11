import numpy as np

def rayTriangleIntersect(point,ray,triangle):
    # point is np array 3x1
    # ray is np array 3x1
    # triangle is np array 3x3
    # return bool, True if intersection
    
    # define an epsilon
    eps = 1e-6
    # break triangle into vertices
    v1, v2, v3 = triangle
    # get edges
    edge1 = v2 - v1
    edge2 = v3 - v1
    # compute
    pvec = np.cross(ray, edge2)
    det = edge1.dot(pvec)
    # check if ray is parallel to triangle
    if abs(det) < eps:
        return False
    # compute ???
    inv_det = 1./ det
    tvec = point - v1
    u = tvec.dot(pvec) * inv_det
    # check ???
    if u < 0. or u > 1.:
        return False
    # compute ???
    qvec = np.cross(tvec, edge1)
    v = ray.dot(qvec) * inv_det
    # check ???
    if v < 0. or u + v > 1.:
        return False
    # compute ???
    t = edge2.dot(qvec) * inv_det
    if t < eps:
        return False
    # if none are tripped then the ray intersects the triangle!
    return True