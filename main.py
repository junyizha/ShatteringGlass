from vpython import *
from numpy import *
from math import *


# Constructing the object Node
class Node:
    def __init__(self, mass, pos, momentum):
        self._mass = mass
        self._pos = pos
        self._momentum = momentum

    def visualize(self):
        sphere(pos=self._pos, radius=10, color=color.red)

    def getMass(self):
        return self._mass

    def getPos(self):
        return self._pos

    def getMomentum(self):
        return self._momentum

    def setPos(self, pos):
        self._pos = pos


# Constructing the object Bond
class Bond:
    def __init__(self, connectedNodes, pos, constants):
        # connectedNodes is a tuple, constants is an array
        self._connectedNodes = connectedNodes
        self._pos = pos
        self._constants = constants

    def visualize(self):
        helix(pos=self._connectedNodes[0].getPos(),
              axis=self._connectedNodes[1].getPos() - self._connectedNodes[0].getPos(),
              radius=0.5, thickness=2, coils=1)


class Floor:
    def __init__(self, altitude):
        self._altitude = altitude

    def visualize(self):
        box(pos=vector(0, self._altitude, 0), length=900, height=1, width=900)


# Constructing the scene
scene2 = canvas(title="Illustration of node", caption="Animated Display", center=vector(0, 0, 0),
                background=color.black)

# Some basic variables
nodes = []
nodesCoordinates = []
bonds = []
bondsCoordinates = []
glassSize = [10, 10, 10]
height = 10

# Creating lattices
for i in range(glassSize[0]):
    for j in range(glassSize[1]):
        for k in range(glassSize[2]):
            nodes.append(Node(10, vector(20 * i, 20 * j + height, 20 * k), 0))
            nodesCoordinates.append([i, j, k])

for i in nodes:
    i.visualize()

# limitshort = 0
# limitlong = 0
# for i in nodesCoordinates:
#     for j in nodesCoordinates:
#         if limitshort == 2 and limitshort == 12:
#             limitshort = 0
#             limitlong = 0
#             break
#         if nodesCoordinates.index(i) < nodesCoordinates.index(j):
#             if (j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 1:
#                 bonds.append(Bond([nodes[nodesCoordinates.index(i)], nodes[nodesCoordinates.index(j)]],
#                                   nodes[nodesCoordinates.index(j)].getPos() - nodes[nodesCoordinates.index(i)].getPos(),
#                                   []))
#                 bondsCoordinates.append([i[0], i[1], i[2], j[0], j[1], j[2]])
#                 limitshort += 1
#             elif (j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 3:
#                 bonds.append(Bond([nodes[nodesCoordinates.index(i)], nodes[nodesCoordinates.index(j)]],
#                                   nodes[nodesCoordinates.index(j)].getPos() - nodes[nodesCoordinates.index(i)].getPos(),
#                                   []))
#                 bondsCoordinates.append([i[0], i[1], i[2], j[0], j[1], j[2]])
#                 limitlong += 1
#
# for i in bonds:
#     i.visualize()

floor1 = Floor(0)

# force analysis
t = 0
dt = 0.1
gravity = 9.8
while t < 100:
    t += dt
    for n in nodes:
        Fgravity = vector(0, -gravity * n.getMass(), 0)
        velocity = Fgravity * dt / n.getMass()
        n.setPos(velocity * dt)
