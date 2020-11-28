from vpython import *
from numpy import *
from math import *


#
# Constructing the object Node
#
class Node:
    def __init__(self, mass, pos, velocity):
        self._mass = mass
        self._pos = pos
        self._velocity = velocity
        self._ball = sphere(pos=self._pos, radius=bondConstants[0] / 2, color=color.red)

    def getMass(self):
        return self._mass

    def getPos(self):
        return self._pos

    def getVelocity(self):
        return self._velocity

    def setPos(self, pos):
        self._pos = pos

    def setVelocity(self, velocity):
        self._velocity = velocity

    def getBall(self):
        return self._ball


#
# Constructing the object Bond
#
class Bond:
    def __init__(self, connectedNodes, pos, constants):
        # connectedNodes is a tuple, constants is an array
        self._connectedNodes = connectedNodes
        self._pos = pos
        self._constants = constants
        self._spring = cylinder(pos=self._connectedNodes[0].getPos(),
                                axis=self._connectedNodes[1].getPos() - self._connectedNodes[0].getPos(),
                                radius=0.0001)

    def getConnectedNodes(self):
        return self._connectedNodes

    def getSpring(self):
        return self._spring

    def setPos(self, pos):
        self._pos = pos

    def getConstants(self):
        return self._constants

    def setConstants(self, constants):
        self._constants = constants


#
# Constructing the object Floor
#
class Floor:
    def __init__(self, altitude):
        self._altitude = altitude
        self._box = box(pos=vector(0, self._altitude, 0), length=9, height=0.001, width=9, color=color.green)


#
# Constructing the scene
#
scene = canvas(title="Shattering Glass", caption="Animated Display", center=vector(0, 0, 0),
               background=color.black)

#
# Some basic variables
#
nodes = []
nodesCoordinates = []
bonds = []
bondsCoordinates = []
glassSize = [3, 3, 3]
nodeMass = 0.0001
height = 0.02
bondConstants = [0.01 / 3, 19, 0, False]  # equilibrium length, spring constant, force, isBroken
bondConstants2 = [0.01 / sqrt(3), 19, 0, False]  # for diagonal bonds

#
# Creating lattices
#
for i in range(glassSize[0]):
    for j in range(glassSize[1]):
        for k in range(glassSize[2]):
            nodes.append(Node(nodeMass, vector(bondConstants[0] * i, bondConstants[0] * j + height,
                                               bondConstants[0] * k), vector(0, 0, 0)))
            nodesCoordinates.append([i, j, k])

#
# Old bond creation
#
limitshort = 0
limitlong = 0
for i in nodesCoordinates:
    for j in nodesCoordinates:
        if limitlong == 4 or limitshort == 12:
            limitshort = 0
            limitlong = 0
            break
        if nodesCoordinates.index(i) < nodesCoordinates.index(j):
            if (j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 1:
                bonds.append(Bond([nodes[nodesCoordinates.index(i)], nodes[nodesCoordinates.index(j)]],
                                  nodes[nodesCoordinates.index(j)].getPos() - nodes[nodesCoordinates.index(i)].getPos(),
                                  bondConstants))
                bondsCoordinates.append([i[0], i[1], i[2], j[0], j[1], j[2]])
                limitshort += 1
            elif (j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 3:
                bonds.append(Bond([nodes[nodesCoordinates.index(i)], nodes[nodesCoordinates.index(j)]],
                                  nodes[nodesCoordinates.index(j)].getPos() - nodes[nodesCoordinates.index(i)].getPos(),
                                  bondConstants2))
                bondsCoordinates.append([i[0], i[1], i[2], j[0], j[1], j[2]])
                limitlong += 1

#
# New bond creation
#
# for i in nodesCoordinates:
#     if ([i[0] + 1, i[1], i[2]] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] + 1, i[1], i[2]] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] + 1, i[1], i[2]])]],
#                           nodes[nodesCoordinates.index([i[0] + 1, i[1], i[2]])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1], i[2]])
#
#     if ([i[0], i[1] + 1, i[2]] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0], i[1] + 1, i[2]] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0], i[1] + 1, i[2]])]],
#                           nodes[nodesCoordinates.index([i[0], i[1] + 1, i[2]])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1] + 1, i[2]])
#
#     if ([i[0], i[1], i[2] + 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0], i[1], i[2] + 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0], i[1], i[2] + 1])]],
#                           nodes[nodesCoordinates.index([i[0], i[1], i[2] + 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1], i[2] + 1])
#
#     if ([i[0] - 1, i[1], i[2]] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] - 1, i[1], i[2]] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] - 1, i[1], i[2]])]],
#                           nodes[nodesCoordinates.index([i[0] - 1, i[1], i[2]])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1], i[2]])
#
#     if ([i[0], i[1] - 1, i[2]] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0], i[1] - 1, i[2]] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0], i[1] - 1, i[2]])]],
#                           nodes[nodesCoordinates.index([i[0], i[1] - 1, i[2]])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1] - 1, i[2]])
#
#     if ([i[0], i[1], i[2] - 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0], i[1], i[2] - 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0], i[1], i[2] - 1])]],
#                           nodes[nodesCoordinates.index([i[0], i[1], i[2] - 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1], i[2] - 1])
#
#     if ([i[0] + 1, i[1] + 1, i[2] + 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] + 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] + 1])]],
#                           nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] + 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] + 1])
#
#     if ([i[0] + 1, i[1] + 1, i[2] - 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] - 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] - 1])]],
#                           nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] - 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] - 1])
#
#     if ([i[0] + 1, i[1] - 1, i[2] + 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] + 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] + 1])]],
#                           nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] + 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] + 1])
#
#     if ([i[0] - 1, i[1] + 1, i[2] + 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] + 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] + 1])]],
#                           nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] + 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] + 1])
#
#     if ([i[0] - 1, i[1] - 1, i[2] + 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] + 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] + 1])]],
#                           nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] + 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] + 1])
#
#     if ([i[0] + 1, i[1] - 1, i[2] - 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] - 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] - 1])]],
#                           nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] - 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] - 1])
#
#     if ([i[0] - 1, i[1] + 1, i[2] - 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] - 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] - 1])]],
#                           nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] - 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] - 1])
#
#     if ([i[0] - 1, i[1] - 1, i[2] - 1] in nodesCoordinates) and (
#             [i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] - 1] not in bondsCoordinates):
#         bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
#                            nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] - 1])]],
#                           nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] - 1])].getPos() - nodes[
#                               nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
#         bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] - 1])


floor = Floor(0)

# test = False
# for i in range(len(bondsCoordinates)):
#     for j in range(i + 1, len(bondsCoordinates)):
#         if i == j:
#             print(i)
#             test = True
#
# print(bondsCoordinates)
# print(test)


#
# force analysis
#
t = 0
dt = 0.0001
gravity = 9.8


while t < 25:
    # sleep(.000001)
    t += dt

    for n in nodes:
        if n.getPos().dot(vector(0, 1, 0)) <= bondConstants[0] / 2:
            n.setVelocity(vector(0, 0, 0) - n.getVelocity())
        n.setPos(n.getPos() + n.getVelocity() * dt)
        n.getBall().pos = n.getPos()  # visualization step

    for b in bonds:
        start = b.getConnectedNodes()[0]
        end = b.getConnectedNodes()[1]
        b.getSpring().axis = end.getPos() - start.getPos()
        b.getSpring().pos = start.getPos()  # visualization step
        force = norm(b.getSpring().axis) * bondConstants[1] * (mag(b.getSpring().axis) - b.getConstants()[0])
        start.setVelocity(start.getVelocity() + force * dt / nodeMass + vector(0, -gravity * nodeMass, 0) * dt / nodeMass -
                          0.02 * start.getVelocity() * dt / nodeMass)
        end.setVelocity(end.getVelocity() - force * dt / nodeMass + vector(0, -gravity * nodeMass, 0) * dt / nodeMass -
                        0.02 * end.getVelocity() * dt / nodeMass)

print("Completed.")
