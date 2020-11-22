from vpython import *
from numpy import *
from math import *


def almostEqual(d1, d2, epsilon=10**-7):
    return (abs(d2 - d1) < epsilon)

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


class Floor:
    def __init__(self, altitude):
        self._altitude = altitude
        self._box = box(pos=vector(0, self._altitude, 0), length=9, height=0.001, width=9, color=color.green)


#
# Constructing the scene
#
scene2 = canvas(title="Illustration of node", caption="Animated Display", center=vector(0, 0, 0),
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
bondConstants = [0.01 / 3, 0.0001, False]  # length, spring constant, isBroken

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
# limitshort = 0
# limitlong = 0
# for i in nodesCoordinates:
#     for j in nodesCoordinates:
#         if limitlong == 4 and limitshort == 12:
#             limitshort = 0
#             limitlong = 0
#             break
#         if nodesCoordinates.index(i) < nodesCoordinates.index(j):
#             if (j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 1:
#                 bonds.append(Bond([nodes[nodesCoordinates.index(i)], nodes[nodesCoordinates.index(j)]],
#                                   nodes[nodesCoordinates.index(j)].getPos() - nodes[nodesCoordinates.index(i)].getPos(),
#                                   bondConstants))
#                 bondsCoordinates.append([i[0], i[1], i[2], j[0], j[1], j[2]])
#                 limitshort += 1
#             elif (j[0] - i[0]) ** 2 + (j[1] - i[1]) ** 2 + (j[2] - i[2]) ** 2 == 3:
#                 bonds.append(Bond([nodes[nodesCoordinates.index(i)], nodes[nodesCoordinates.index(j)]],
#                                   nodes[nodesCoordinates.index(j)].getPos() - nodes[nodesCoordinates.index(i)].getPos(),
#                                   bondConstants))
#                 bondsCoordinates.append([i[0], i[1], i[2], j[0], j[1], j[2]])
#                 limitlong += 1

#
# New bond creation
#
for i in nodesCoordinates:
    if ([i[0] + 1, i[1], i[2]] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1], i[2]] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1], i[2]])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1], i[2]])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1], i[2]])

    if ([i[0], i[1] + 1, i[2]] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0], i[1] + 1, i[2]] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0], i[1] + 1, i[2]])]],
                          nodes[nodesCoordinates.index([i[0], i[1] + 1, i[2]])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1] + 1, i[2]])

    if ([i[0], i[1], i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0], i[1], i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0], i[1], i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0], i[1], i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1], i[2] + 1])

    if ([i[0] - 1, i[1], i[2]] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1], i[2]] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1], i[2]])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1], i[2]])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1], i[2]])

    if ([i[0], i[1] - 1, i[2]] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0], i[1] - 1, i[2]] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0], i[1] - 1, i[2]])]],
                          nodes[nodesCoordinates.index([i[0], i[1] - 1, i[2]])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1] - 1, i[2]])

    if ([i[0], i[1], i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0], i[1], i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0], i[1], i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0], i[1], i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0], i[1], i[2] - 1])

    if ([i[0] + 1, i[1] + 1, i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] + 1])

    if ([i[0] + 1, i[1] + 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] - 1])

    if ([i[0] + 1, i[1] - 1, i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] + 1])

    if ([i[0] - 1, i[1] + 1, i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] + 1])

    if ([i[0] - 1, i[1] - 1, i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] + 1])

    if ([i[0] + 1, i[1] - 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] - 1])

    if ([i[0] - 1, i[1] + 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] - 1])

    if ([i[0] - 1, i[1] - 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] - 1])

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
dt = 0.001
gravity = 9.8
Fnet = vector(0, -gravity * nodeMass, 0)
acceleration = Fnet / nodeMass

while t < 25:
    sleep(0.1)
    t += dt
    for n in nodes:
        X = nodesCoordinates[nodes.index(n)][0]
        Y = nodesCoordinates[nodes.index(n)][1]
        Z = nodesCoordinates[nodes.index(n)][2]
        tempnodes = [[X, Y, Z, X + 1, Y, Z],
                     [X, Y, Z, X, Y + 1, Z], [X, Y, Z, X, Y, Z + 1], [X, Y, Z, X - 1, Y, Z],
                     [X, Y, Z, X, Y - 1, Z], [X, Y, Z, X, Y, Z - 1], [X, Y, Z, X + 1, Y + 1, Z + 1],
                     [X, Y, Z, X + 1, Y - 1, Z + 1],
                     [X, Y, Z, X + 1, Y + 1, Z - 1], [X, Y, Z, X - 1, Y + 1, Z + 1], [X, Y, Z, X - 1, Y - 1, Z + 1],
                     [X, Y, Z, X - 1, Y + 1, Z - 1], [X, Y, Z, X - 1, Y - 1, Z - 1], [X, Y, Z, X + 1, Y - 1, Z - 1]]

        if n.getPos().dot(vector(0, 1, 0)) <= bondConstants[0]/2:
            n.setVelocity(-n.getVelocity() * 0.8)
        # for m in nodes:
        #     if abs(mag(n.getPos()) - mag(m.getPos())) >= 0.1 and mag(n.getPos())**2 + mag(m.getPos())**2 < 99.9:
        #         n.setMomentum(-n.getMomentum())
        #         m.setMomentum(-m.getMomentum())
        # for j in tempnodes:
        #     if j in bondsCoordinates:
        #         #if (abs((j[0] - j[3])) + abs((j[1] - j[4])) + abs((j[2] - j[5]))) == 1:
        #         if almostEqual((abs((j[0] - j[3])) + abs((j[1] - j[4])) + abs((j[2] - j[5]))), 1, epsilon=10**-7):
        #             cNodes = bonds[bondsCoordinates.index(j)].getConnectedNodes()
        #             springLength = n.getPos() - cNodes[1].getPos()
        #             Fnet += (bondConstants[0] * norm(springLength) - springLength) * bondConstants[1]
                # else:
                #     cNodes = bonds[bondsCoordinates.index(j)].getConnectedNodes()
                #     springLength = n.getPos() - cNodes[1].getPos()
                #     Fnet += (sqrt(3) * bondConstants[0] * norm(springLength) - springLength) * bondConstants[1]

        if n.getPos().dot(vector(0, 1, 0)) <= bondConstants[0] / 2:
            Fnet += norm(vector(0, 1, 0)) * (bondConstants[0] - n.getPos().dot(vector(0, 1, 0))) * 9999
        n.setPos(n.getPos() + n.getVelocity() * dt + 0.5 * acceleration * dt ** 2)  # update position
        n.getBall().pos = n.getPos()  # visualization step

        newAcceleration = Fnet / nodeMass
        n.setVelocity(n.getVelocity() + 0.5 * (acceleration + newAcceleration) * dt)  # update velocity

        acceleration = newAcceleration
        Fnet = vector(0, -gravity * nodeMass, 0)

        for j in tempnodes:
            if j in bondsCoordinates:
                bonds[bondsCoordinates.index(j)].getSpring().axis = nodes[nodesCoordinates.index(
                    [j[3], j[4], j[5]])].getPos() - n.getPos()
                bonds[bondsCoordinates.index(j)].getSpring().pos = n.getPos()
