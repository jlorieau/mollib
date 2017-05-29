import timeit

no_points = 2000

setup = """
import numpy as np
from collections import namedtuple
from operator import itemgetter
from math import sqrt
from pprint import pformat, pprint

np.random.seed(0)
point_list = [np.random.rand(3)*20. -10. for i in range({no_points})]

class Node(namedtuple('Node', 'location left_child right_child')):
    def __repr__(self):
        return pformat(tuple(self))

def kdtree(point_list, depth=0):
    try:
        k = len(point_list[0]) # assumes all points have the same dimension
    except IndexError as e: # if not point_list:
        return None
    # Select axis based on depth so that axis cycles through all valid values
    axis = depth % k

    # Sort point list and choose median as pivot element
    point_list.sort(key=itemgetter(axis))
    median = len(point_list) // 2 # choose median

    # Create node and construct subtrees
    return Node(
        location=point_list[median],
        left_child=kdtree(point_list[:median], depth + 1),
        right_child=kdtree(point_list[median + 1:], depth + 1)
    )

def dist(p1, p2):
    return sqrt(sum([(i-j)**2 for i,j in zip(p1, p2)]))

def find_within(point, closest_list, kdtree, cutoff= 1.0, depth=0):
    "returns atoms within cutoff and the remaining kdtree nodes to search"

    dist0 =dist(point, kdtree[0])
    if dist0 < cutoff:
        closest_list.append(kdtree[0])

    axis = depth % len(point)

    pt = point[axis]
    pt = kdtree[0][axis]

    l1 = pt + cutoff*1.
    l2 = pt - cutoff*1.


    if kdtree[1] is not None:
        find_within(point, closest_list, kdtree[1], cutoff, depth+1)
    if kdtree[2] is not None:
        find_within(point, closest_list, kdtree[2], cutoff, depth+1)

point = (-3.2, 4.5, -1.5)
closest = []
tree = kdtree(point_list)
""".format(no_points=no_points)

setup2 = """
from scipy import spatial
import numpy as np
np.random.seed(0)
points = [np.random.rand(3)*20. -10. for i in range({no_points})]

tree = spatial.cKDTree(points)
""".format(no_points=no_points)

setup3 = """
import numpy as np
from collections import namedtuple
from operator import itemgetter
from math import sqrt
from pprint import pformat, pprint

np.random.seed(0)
point_list = [np.random.rand(3)*20. -10. for i in range({no_points})]

class Node(namedtuple('Node', 'location left_child right_child')):
    def __repr__(self):
        return pformat(tuple(self))

def kdtree(point_list, depth=0):
    try:
        k = len(point_list[0]) # assumes all points have the same dimension
    except IndexError as e: # if not point_list:
        return None
    # Select axis based on depth so that axis cycles through all valid values
    axis = depth % k

    # Sort point list and choose median as pivot element
    point_list.sort(key=itemgetter(axis))
    median = len(point_list) // 2 # choose median

    # Create node and construct subtrees
    return Node(
        location=point_list[median],
        left_child=kdtree(point_list[:median], depth + 1),
        right_child=kdtree(point_list[median + 1:], depth + 1)
    )

def dist(p1, p2):
    return sqrt(sum([(i-j)**2 for i,j in zip(p1, p2)]))

def find_within(point, closest_list, kdtree, cutoff= 1.0, depth=0):
    "returns atoms within cutoff and the remaining kdtree nodes to search"

    dist0 =dist(point, kdtree[0])
    if dist0 < cutoff:
        closest_list.append(kdtree[0])

    axis = depth % len(point)

    #pt = kdtree[0][axis]
    pt = point[axis]

    l1 = pt + cutoff*1.
    l2 = pt - cutoff*1.

    go_left = False
    go_right = False

    if kdtree[1] is not None:
        if kdtree[1][0][axis] < l1:
            go_left |= True
        if kdtree[1][0][axis] > l2:
            go_left |= True
    if kdtree[2] is not None:
        if kdtree[2][0][axis] < l1:
            go_right |= True
        if kdtree[2][0][axis] > l2:
            go_right |= True

    if go_left:
        find_within(point, closest_list, kdtree[1], cutoff, depth+1)
    if go_right:
        find_within(point, closest_list, kdtree[2], cutoff, depth+1)
point = (-3.2, 4.5, -1.5)
closest = []
tree = kdtree(point_list)
""".format(no_points=no_points)

time = timeit.timeit('find_within(point, closest, tree, cutoff=2.0)',
                     setup=setup,
                     number=no_points)
print ('kd-tree: {} ms/iter'.format(time * 1000. / no_points))

time = timeit.timeit('find_within(point, closest, tree, cutoff=2.0)',
                     setup=setup3,
                     number=no_points)
print ('kd-tree non-lin search: {} ms/iter'.format(time * 1000. / no_points))

time = timeit.timeit('[i for i in point_list if dist(i, point) < 2.0]',
                     setup=setup,
                     number=no_points)
print ('linear: {} ms/iter'.format(time * 1000. / no_points))

time = timeit.timeit(setup2, '', number=1)
print ('scipy/setup: {} ms/iter'.format(time))

time = timeit.timeit('tree.query_ball_point((-3.2, 4.5, -1.5), 2.0)', setup2,
                     number=no_points)
print ('scipy/kdtree: {} ms/iter'.format(time * 1000. / no_points))
