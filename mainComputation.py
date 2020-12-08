import csv

from vpython import *
from numpy import *
from math import *


#
# Constructing the object Sphere
#
class Sphere:
    def __init__(self, pos):
        self._pos = pos

    def getPos(self):
        return self._pos

    def setPos(self, newpos):
        self._pos = newpos


#
# Constructing the object Cylinder
#
class Cylinder:
    def __init__(self, pos, axis):
        self._pos = pos
        self._axis = axis

    def getAxis(self):
        return self._axis

    def setAxis(self, newAxis):
        self._axis = newAxis

    def getPos(self):
        return self._pos

    def setPos(self, updatePos):
        self._pos = updatePos


#
# Constructing the object newBox
#
class newBox:
    def __init__(self, pos, up):
        self._pos = pos
        self._up = up

    def getPos(self):
        return self._pos

    def getUp(self):
        return self._up

    def setPos(self, pos):
        self._pos = pos

    def setUp(self, up):
        self._up = up


#
# Constructing the object Node
#
class Node:
    def __init__(self, mass, pos, velocity):
        self._mass = mass
        self._pos = pos
        self._velocity = velocity
        self._ball = Sphere(self._pos)

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
        self._spring = Cylinder(self._connectedNodes[0].getPos(),self._connectedNodes[1].getPos() - self._connectedNodes[0].getPos())
        self._lastforce = vector(0, 0, 0)

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

    def getlastForce(self):
        return self._lastforce

    def setlastForce(self, inputforce):
        self._lastforce = inputforce


#
# Constructing the object Floor
#
class Floor:
    def __init__(self, altitude):
        self._altitude = altitude
        self._box = newBox(vector(0, self._altitude, 0), norm(vector(0.1, 1, 0)))

    def shortestDistance(self, posvector):
        d = abs((self._box.getUp().x * posvector.x + self._box.getUp().y * posvector.y + self._box.getUp().z * posvector.z))
        e = (math.sqrt(
            self._box.getUp().x * self._box.getUp().x + self._box.getUp().y * self._box.getUp().y + self._box.getUp().z * self._box.getUp().z))
        return d / e

    def getNormalVector(self):
        return self._box.getUp()


#
# Some basic variables
#
nodes = []
nodesCoordinates = []
bonds = []
bondsCoordinates = []
glassSize = [3, 3, 3]
nodeMass = 0.0001
height = 0.05
bondConstants = [0.01 / 3, 999, 0, False]  # equilibrium length, spring constant, force, isBroken
bondConstants2 = [0.01 / sqrt(3), 999 * sqrt(3), 0, False]  # for diagonal bonds


#
# Creating lattices
#
for i in range(glassSize[0]):
    for j in range(glassSize[1]):
        for k in range(glassSize[2]):
            nodes.append(Node(nodeMass, vector(bondConstants[0] * i, bondConstants[0] * j + height,
                                               bondConstants[0] * k), vector(0, 0, 0)), )
            nodesCoordinates.append([i, j, k])


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
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] + 1])

    if ([i[0] + 1, i[1] + 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1] + 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] + 1, i[2] - 1])

    if ([i[0] + 1, i[1] - 1, i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] + 1])

    if ([i[0] - 1, i[1] + 1, i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] + 1])

    if ([i[0] - 1, i[1] - 1, i[2] + 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] + 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] + 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] + 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] + 1])

    if ([i[0] + 1, i[1] - 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] + 1, i[1] - 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] + 1, i[1] - 1, i[2] - 1])

    if ([i[0] - 1, i[1] + 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] + 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] + 1, i[2] - 1])

    if ([i[0] - 1, i[1] - 1, i[2] - 1] in nodesCoordinates) and (
            [i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] - 1] not in bondsCoordinates):
        bonds.append(Bond([nodes[nodesCoordinates.index([i[0], i[1], i[2]])],
                           nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] - 1])]],
                          nodes[nodesCoordinates.index([i[0] - 1, i[1] - 1, i[2] - 1])].getPos() - nodes[
                              nodesCoordinates.index([i[0], i[1], i[2]])].getPos(), bondConstants2))
        bondsCoordinates.append([i[0], i[1], i[2], i[0] - 1, i[1] - 1, i[2] - 1])

floor = Floor(0)

#
# force analysis
#
t = 0
dt = 0.00001
gravity = 9.8
count = 0
isCollided = []
oldAcceleration = vector(0, -gravity, 0)

for i in range(len(nodes)):
    isCollided.append(False)

#
# Main While Loop
#
f = open("output.csv", "w+")
row = []
csvwriter = csv.writer(f)
heading = []
for i in range(glassSize[0] * glassSize[1] * glassSize[2]):
    tempStr = ("Node " + str(i))
    heading.append(tempStr)
csvwriter.writerow(heading)
count = 100

while t < 25:
    # sleep(.0000005)
    t += dt

    #
    # Mechanics of nodes
    #
    for n in nodes:
        nCoordinates = nodesCoordinates[nodes.index(n)]

        #
        # Colliding with the floor
        #
        if floor.shortestDistance(n.getPos()) <= bondConstants[0] / 2:
            projectedY = n.getVelocity().dot(floor.getNormalVector()) / mag(floor.getNormalVector()) * norm(
                floor.getNormalVector())
            projectedX = n.getVelocity() - projectedY
            newProjectedY = -projectedY
            n.setVelocity(1.0 * (newProjectedY + projectedX))

        #
        # Collision Detection
        #
        for m in nodes:
            mCoordinates = nodesCoordinates[nodes.index(m)]
            connected = False
            if [nCoordinates[0], nCoordinates[1], nCoordinates[2], mCoordinates[0], mCoordinates[1],
                mCoordinates[2]] in bondsCoordinates:
                mnCoordinates = [nCoordinates[0], nCoordinates[1], nCoordinates[2], mCoordinates[0], mCoordinates[1],
                                 mCoordinates[2]]
                connected = True
            elif [mCoordinates[0], mCoordinates[1], mCoordinates[2], nCoordinates[0], nCoordinates[1],
                  nCoordinates[2]] in bondsCoordinates:
                mnCoordinates = [mCoordinates[0], mCoordinates[1], mCoordinates[2], nCoordinates[0], nCoordinates[1],
                                 nCoordinates[2]]
                connected = True
            if (not connected or connected and bonds[bondsCoordinates.index(mnCoordinates)].getConstants()[3]) and \
                    0.000001 <= mag(n.getPos() - m.getPos()) < bondConstants[0] and not isCollided[nodes.index(m)]:
                isCollided[nodes.index(m)] = True
                distance = mag(m.getPos() - n.getPos())
                vCenterN = (m.getPos() - n.getPos()).dot(n.getVelocity()) / distance * norm(m.getPos() - n.getPos())
                vCenterM = (n.getPos() - m.getPos()).dot(m.getVelocity()) / distance * norm(n.getPos() - m.getPos())
                n.setVelocity(n.getVelocity() - vCenterN + vCenterM)
                m.setVelocity(m.getVelocity() - vCenterM + vCenterN)
                if (m.getPos() - n.getPos()).dot(n.getVelocity()) / mag(m.getPos() - n.getPos()) < 0:
                    isCollided[nodes.index(m)] = False

        n.setVelocity(n.getVelocity() + 0.5 * dt * oldAcceleration)
        n.setPos(n.getPos() + dt * n.getVelocity())
        n.setVelocity(n.getVelocity() + 0.5 * dt * oldAcceleration)
        # n.setVelocity(n.getVelocity() + vector(0, -gravity * nodeMass, 0) * dt / nodeMass)
        # n.setPos(n.getPos() + n.getVelocity() * dt)
        # n.getBall().pos = n.getPos()  # visualization step

        row.append(str(n.getPos()))

    # count += 1

    #
    # Converting to csv
    #
    # if count == 100:
    #     csvwriter.writerow(row)
    #     count = 0
    csvwriter.writerow(row)
    row = []

    #
    # Mechanics of bonds
    #
    for b in bonds:
        if b.getConstants()[3]:
            continue

        start = b.getConnectedNodes()[0]
        end = b.getConnectedNodes()[1]
        b.getSpring().setAxis(end.getPos() - start.getPos())
        b.getSpring().setPos(start.getPos())
        springStretch = mag(b.getSpring().getAxis()) - b.getConstants()[0]
        force = norm(b.getSpring().getAxis()) * b.getConstants()[1] * springStretch

        # zeta = 2 * sqrt(nodeMass / b.getConstants()[1])
        # start.setVelocity(start.getVelocity() + force * dt / nodeMass - 0.5 * zeta * norm(start.getVelocity()) * mag(
        #     start.getVelocity() - end.getVelocity()) * dt / nodeMass)
        # end.setVelocity(end.getVelocity() - force * dt / nodeMass - 0.5 * zeta * norm(end.getVelocity()) * mag(
        #     end.getVelocity() - start.getVelocity()) * dt / nodeMass)

        if t == dt:
            b.setlastForce(force)
        zeta = 2 * sqrt(nodeMass / b.getConstants()[1])
        start.setVelocity(start.getVelocity() + 0.5 * dt * (
                b.getlastForce() - 0.5 * zeta * norm(start.getVelocity()) * mag(
            start.getVelocity() - end.getVelocity()) / nodeMass))
        start.setVelocity(start.getVelocity() + 0.5 * dt * (force - 0.5 * zeta * norm(start.getVelocity()) * mag(
            start.getVelocity() - end.getVelocity())) / nodeMass)
        end.setVelocity(end.getVelocity() - 0.5 * dt * (b.getlastForce() + 0.5 * zeta * norm(end.getVelocity()) * mag(
            end.getVelocity() - start.getVelocity())) / nodeMass)
        end.setVelocity(end.getVelocity() - 0.5 * dt * (force + 0.5 * zeta * norm(end.getVelocity()) * mag(
            end.getVelocity() - start.getVelocity()) / nodeMass))
        b.setlastForce(force)

        #
        # Breakage modelling
        #
        breakageThreshold = 0.05
        if springStretch >= breakageThreshold * b.getConstants()[0]:
            b.setConstants([b.getConstants()[0], b.getConstants()[1], b.getConstants()[2], True])

f.close()
print("Computation Completed.")
