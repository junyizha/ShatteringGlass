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


# Constructing the scene
scene2 = canvas(title="Illustration of node", caption="Animated Display", center=vector(0, 0, 0),
                background=color.black)

# Some basic arrays
nodes = []
nodesCoordinates = []
bonds = []
bondsCoordinates = []

# Creating lattices
for i in range(3):
    for j in range(3):
        for k in range(3):
            nodes.append(Node(10, vector(30 * i, 30 * j, 30 * k), 0))
            nodesCoordinates.append([i, j, k])

for i in nodes:
    i.visualize()

for i in nodesCoordinates:
    for j in nodesCoordinates:
        if nodesCoordinates.index(i) < nodesCoordinates.index(j) \
                and ((j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 1
                or (j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 3):
            bonds.append(Bond([nodes[nodesCoordinates.index(i)], nodes[nodesCoordinates.index(j)]],
                              nodes[nodesCoordinates.index(j)].getPos() - nodes[nodesCoordinates.index(i)].getPos(),
                              []))
            bondsCoordinates.append([i[0], i[1], i[2], j[0], j[1], j[2]])

for i in bonds:
    i.visualize()

print(nodesCoordinates)
print(bondsCoordinates)
