import numpy as np

def random(d,n):
    # d is the dimension
    # n is the number of samples
    return np.random.rand(n,d)

def poissonDisk(d,n):
    rmin = 1/(n**(1/d)) # reasonable r_min to fill the domain
    pts = np.empty((n,d)) # Initialize the array, we will trim later
    pts[0,:] = np.random.rand(d) # Place the first point
    for i in range(1,n):
        dist = 0
        attempts = 0
        while (dist < rmin):
            xnew = np.random.rand(d)
            dist = np.min(np.linalg.norm(xnew - pts[:i,:], axis = 1))
            if (attempts > 1000):
                return pts[:i,:]
            attempts += 1
        pts[i,:] = xnew
    return pts

def poissonDiskDecay(d,n):
    rmin = 1/(n**(1/d)) # reasonable r_min to fill the domain
    pts = np.empty((n,d)) # Initialize the array, we will trim later
    pts[0,:] = np.random.rand(d) # Place the first point
    for i in range(1,n):
        dist = 0
        attempts = 0
        while (dist < rmin):
            xnew = np.random.rand(d)
            dist = np.min(np.linalg.norm(xnew - pts[:i,:], axis = 1))
            if (attempts > 1000):
                rmin *= 0.95
            attempts += 1
        pts[i,:] = xnew
    return pts

def halton(d,n):
    pts = np.empty((n,d))
    if d == 1:
        for i, x in zip(range(n), haltonGenerator(2)):
            pts[i,0] = x
    elif d == 2:
        for i, x, y in zip(range(n), haltonGenerator(2), haltonGenerator(3)):
            pts[i,0] = x
            pts[i,1] = y
    elif d == 3:
        for i, x, y, z in zip(range(n), haltonGenerator(2), haltonGenerator(3), haltonGenerator(5)):
            pts[i,0] = x
            pts[i,1] = y
            pts[i,2] = z
    else: # d > 3
        raise TypeError('Invalid Dimension (d > 3)')
    return pts

def haltonGenerator(b):
    n, d = 0, 1
    while True:
        x = d - n
        if x == 1:
            n = 1
            d *= b
        else:
            y = d // b
            while x <= y:
                y //= b
            n = (b + 1) * y - x
        yield n / d