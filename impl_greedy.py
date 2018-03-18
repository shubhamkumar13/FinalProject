from math import sqrt
import numpy as np
from scipy.spatial import ConvexHull
from bokeh.plotting import figure, output_file, show
from random import randint
from matplotlib import collections as mc
import collections as cl


class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return "( %s, %s)" % (self.x, self.y)


class Line(object):

    def __init__(self, p1, p2):
        self.point1 = p1
        self.point2 = p2

    def Length(self):
        l = (self.point1.x - self.point2.x)**2 + \
            (self.point1.y - self.point2.y)**2
        l = sqrt(l)
        return l

    def Slope(self):
        numerator = self.point2.y - self.point1.y
        denominator = self.point2.x - self.point1.x
        if denominator == 0:
            m = float("Inf")
        else:
            m = numerator / denominator
        return m

    def Y_intercept(self):
        m = self.Slope()
        return (self.point1.y - m*self.point1.x)

    def X_intercept(self):
        return self.point1.x

    def __repr__(self):
        m = self.Slope()
        c = self.Y_intercept()
        if m == float("inf"):
            c = self.X_intercept()
            return "( x = %s )\n" % (c)
        elif m is 0:
            return "( y = %s )\n" % (c)
        else:
            return "( y = %sx + %s )\n" % (m, c)


class Circle(object):
    """docstring for Circle"""

    def __init__(self, p0, p=None, radius=None):
        self.point_origin = p0
        self.point_end = p
        self.radius = radius

    def Radius(self):
        if radius is None:
            r = (self.point_origin.x - self.point_end.x)**2 + \
                (self.point_origin.y - self.point_end.y)**2
            r = sqrt(r)
            return r
        else:
            return self.radius

    def __repr__(self):
        return "( x - %s )^2 + ( y - %s )^2 = %s^2\n" % (self.point_origin.x, self.point_origin.y, radius)



def distPoint_Point(point1, point2):
    return sqrt(square(point1.x - point2.x) + square(point1.y - point2.y))


def square(x):
    return x**2


def distPoint_Line(point: Point, line: Line):
    '''
    distance between point and line
    '''
    if line.Length() == 0:
        return distPoint_Point(point, line.point1)
    else:
        t = ((point.x - line.point1.x)*(line.point2.x - line.point1.x) +
             (point.y - line.point1.y)*(line.point2.y - line.point1.y)) / line.Length()
        t = max(0, min(1, t))
        new_point = Point(line.point1.x + t*(line.point2.x - line.point1.x),
                          line.point1.y + t*(line.point2.y - line.point1.y))
        return distPoint_Point(point, new_point)


def isIntersection(circle1: Circle, circle2: Circle = None, line: Line = None):
    '''
        checks if a line segment is inside
        by using calculating 
        --if line segment points inside the circle
        --if line segment distance from center < radius

        if above 2 satisfied
            line segment inside the circle
        
        create a dict for (index of circle) -> [all(index of line)]
    '''
    if line is None and circle2 is not None:
        line = Line(circle1.point_origin, circle2.point_origin)
        if line.Length() < circle1.Radius() + circle2.Radius():
            return True
        else:
            return False
    elif line is not None and circle2 is None:
        dist_centre_line = distPoint_Line(circle1.point_origin, line)
        if dist_centre_line < circle1.Radius():
            return True
        else:
            return False
    else:
        return False


'''
def convert_list(points, index):
    x = points[index]
    x = np.asarray(x)
    print(x)
    return x
'''
'''
def convert_ndarray(x, y):
    print(x)
    print(y)
'''


def convert_points_matrix(points):
    '''
    convert the points object into 
    a 2d numpy array
    '''
    temp_list = list()
    matrix = list()
    for p in points:
        temp_list.append(p.x)
        temp_list.append(p.y)
        matrix.append(temp_list)
        temp_list = []

    matrix = np.array(matrix)
    return matrix

class lineCollection(object):
    def __init__(self, lines):
        self.lines = lines
    
    def NumLines(self):
        c = len(self.lines)
        return c

    def DictLines(self):
        c = self.NumLines()
        key = [i for i in range(c)]
        value = [line for line in self.lines]
        dict_lines = cl.OrderedDict(zip(key, value))
        return dict_lines

class circleCollection(object):
    def __init__(self, circles):
        self.circles = circles
    
    def NumCircles(self):
        c = len(self.circles)
        return c
    
    def DictCircles(self):
        pass



def Random_Lines(points):
    '''
        input points in pairs which form lines
        since points are genrated randomly
        need to take point pairs
        that would serve as input lines
    '''
    Lines = list()
    for index in range(1, len(points)):
        temp_lines = Line(points[index-1], points[index])
        Lines.append(temp_lines)
    return Lines


def Random_Circles(points, radius=5):
    '''
        input points not pairs of points
        this means generate points
        then the unit radii would take care
    '''
    circles = list()
    for p in points:
        circles.append(Circle(p0=p, radius=radius))

    return circles


def display_encompass_points(points):
    """
        As stated we need to surround 
        the total lines generated
        so we need a 2d plane of sorts
        this means we need to apply
        Convex Hull
        And we can get back the vertices
        which are part of the convex hull
    """
    hull = ConvexHull(points)

    import matplotlib.pyplot as plt
    plt.plot(points[:, 0], points[:, 1], '+')
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

    plt.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'r--', lw=2)
    plt.plot(points[hull.vertices[0], 0], points[hull.vertices[0], 1])
    plt.show()


def display_lines(lines):
    '''
        create a display of lines
    '''


def display_circles(circles):
    '''
        create a display of circles
    '''
    pass
"""
def LengthIntersection(line:Line, circle:Circle):
    '''
        calculate the length of the line
        intersected by the circle
    '''
    dist_point_line = distPoint_Line(circle.point_origin, line)
    radius = circle.Radius()
    segment = radius*radius - dist_point_line*dist_point_line
    segment = sqrt(segment)
    segment = 2*segment
    return segment

def SumLengthIntersection(lines, circle):
    '''
        calculate delta which sum of all lines
        which intersect with one of the circle
    '''
    delta = 0
    for line in lines:
        delta = delta + LengthIntersection(line, circle)
    return delta

def GetDelta(lines, circles):
    '''
        mapping for circle -> delta value
        return the sorted order
        in decreasing order
    '''
    dict_circle_delta = dict()
    for circle in circles:
        delta = SumLengthIntersection(lines, circle)
        dict_circle_delta[circle] = delta
    
    dict_circle_delta = cl.OrderedDict(dict_circle_delta)
    dict_circle_delta = cl.OrderedDict(sorted(dict_circle_delta.items(), key=lambda t: t[1]))
    return dict_circle_delta

def ChangeDeltaDict(delta_dict, delta_value):
    for key,value in delta_dict.items():
        pass
    pass
"""

def disk_cover_algorithm(circles, lines):
    pass


# Random Points generated using poisson's distribution
def lonely(p, X, r):
    m = X.shape[1]
    x0, y0 = p
    x = y = np.arange(-r, r)
    x = x + x0
    y = y + y0

    u, v = np.meshgrid(x, y)

    u[u < 0] = 0
    u[u >= m] = m-1
    v[v < 0] = 0
    v[v >= m] = m-1

    return not np.any(X[u[:], v[:]] > 0)


def Random_Points(m=2500, r=200, k=30):
    '''
        this creates random points, 
        the random points are returned 
        as an np array
    '''
    # m = extent of sample domain
    # r = minimum distance between points
    # k = samples before rejection
    active_list = []

    # step 0 - initialize n-d background grid
    X = np.ones((m, m))*-1

    # step 1 - select initial sample
    x0, y0 = np.random.randint(0, m), np.random.randint(0, m)
    active_list.append((x0, y0))
    X[active_list[0]] = 1

    # step 2 - iterate over active list
    while active_list:
        i = np.random.randint(0, len(active_list))
        rad = np.random.rand(k)*r+r
        theta = np.random.rand(k)*2*np.pi

        # get a list of random candidates within [r,2r] from the active point
        candidates = np.round(
            (rad*np.cos(theta)+active_list[i][0], rad*np.sin(theta)+active_list[i][1])).astype(np.int32).T

        # trim the list based on boundaries of the array
        candidates = [(x, y) for x, y in candidates if x >=
                      0 and y >= 0 and x < m and y < m]

        for p in candidates:
            if X[p] < 0 and lonely(p, X, r):
                X[p] = 1
                active_list.append(p)
                break
        else:
            del active_list[i]

    '''
    Converts N x 2 array
    to
    A point class
    '''
    X = np.where(X > 0)
    Points = list()
    rangeX = X[0]
    rangeY = X[1]

    for x, y in zip(rangeX, rangeY):
        temp_point = Point(x, y)
        Points.append(temp_point)

    return Points


Points_Lines = Random_Points(250, 10, 10)
Lines = Random_Lines(Points_Lines)
Points_Circles = Random_Points(250, 10, 10)
radius = 5
Circles = Random_Circles(Points_Circles, radius)
