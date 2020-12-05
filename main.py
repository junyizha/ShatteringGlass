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
        self._box = box(pos=vector(0, self._altitude, 0), up=norm(vector(0.1, 1, 0)), length=0.4, height=0.001,
                        width=0.1, color=color.green)

    def shortestDistance(self, posvector):
        d = abs((self._box.up.x * posvector.x + self._box.up.y * posvector.y + self._box.up.z * posvector.z))
        e = (math.sqrt(self._box.up.x * self._box.up.x + self._box.up.y * self._box.up.y + self._box.up.z * self._box.up.z))
        return d / e

    def getNormalVector(self):
        return self._box.up


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
glassSize = [2, 2, 2]
nodeMass = 0.0001
height = 0.05
bondConstants = [0.01 / 3, 9, 0, False]  # equilibrium length, spring constant, force, isBroken
bondConstants2 = [0.01 / sqrt(3), 9 * sqrt(3), 0, False]  # for diagonal bonds


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
#         # if limitlong == 4 or limitshort == 12:
#         #     limitshort = 0
#         #     limitlong = 0
#         #     break
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
#                                   bondConstants2))
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
dt = 0.00001
gravity = 9.8
count = 0
isCollided = []

for i in range(len(nodes)):
    isCollided.append(False)

#
# Main While Loop
#
while t < 25:
    # sleep(.0000005)
    t += dt

    #
    # Keyboard events
    #
    # k = keysdown()
    # if 'w' in k:
    #     scene.camera.pos = scene.camera.pos + vector(0, .001, 0)
    #     scene.forward = scene.center - scene.camera.pos
    # if 's' in k:
    #     scene.camera.pos = scene.camera.pos - vector(0, .001, 0)
    #     scene.forward = scene.center - scene.camera.pos
    # if 'd' in k:
    #     scene.camera.pos = scene.camera.pos + vector(.001, 0, 0)
    #     scene.forward = scene.center - scene.camera.pos
    # if 'a' in k:
    #     scene.camera.pos = scene.camera.pos - vector(.001, 0, 0)
    #     scene.forward = scene.center - scene.camera.pos
    # if 'h' in k:
    #     scene.camera.pos = scene.camera.pos + vector(0, 0, .001)
    #     scene.forward = scene.center - scene.camera.pos
    # if 'n' in k:
    #     scene.camera.pos = scene.camera.pos - vector(0, 0, .001)
    #     scene.forward = scene.center - scene.camera.pos

    #
    # Mechanics of nodes
    #
    for n in nodes:

        #
        # Colliding with the floor
        #
        if floor.shortestDistance(n.getPos()) <= bondConstants[0] / 2:
            # print(floor.shortestDistance(n.getPos()))
            projectedY = n.getVelocity().dot(floor.getNormalVector()) / mag(floor.getNormalVector()) * norm(
                floor.getNormalVector())
            projectedX = n.getVelocity() - projectedY
            newProjectedY = -projectedY
            n.setVelocity(1.0 * (newProjectedY + projectedX))

        #
        # if n.getPos().y <= bondConstants[0] / 2:
        #     n.setVelocity(n.getVelocity() + (
        #                 vector(0, -gravity * nodeMass, 0) + vector(0, -(n.getPos().y - bondConstants[0] / 2) * 99999,
        #                                                            0)) * dt / nodeMass)
        # else:
        #     n.setVelocity(n.getVelocity() + vector(0, -gravity * nodeMass, 0) * dt / nodeMass)

        # n.setVelocity(vector(0, 0, 0) - n.getVelocity())

        #
        # Collision Detection
        #
        # for m in nodes:
        #     if 0.000001 <= mag(n.getPos() - m.getPos()) < bondConstants[0] and not isCollided[nodes.index(m)]:
        #         isCollided[nodes.index(m)] = True
        #         distance = mag(m.getPos() - n.getPos())
        #         vCenterN = (m.getPos() - n.getPos()).dot(n.getVelocity()) / distance * norm(m.getPos() - n.getPos())
        #         vCenterM = (n.getPos() - m.getPos()).dot(m.getVelocity()) / distance * norm(n.getPos() - m.getPos())
        #         n.setVelocity(n.getVelocity() - vCenterN + vCenterM)
        #         m.setVelocity(m.getVelocity() - vCenterM + vCenterN)
        #         print("Collided")
        #         if (m.getPos() - n.getPos()).dot(n.getVelocity()) / mag(m.getPos() - n.getPos()) < 0:
        #             isCollided[nodes.index(m)] = False

        n.setVelocity(n.getVelocity() + vector(0, -gravity * nodeMass, 0) * dt / nodeMass)
        n.setPos(n.getPos() + n.getVelocity() * dt)
        n.getBall().pos = n.getPos()  # visualization step

    #
    # Mechanics of bonds
    #
    for b in bonds:
        start = b.getConnectedNodes()[0]
        end = b.getConnectedNodes()[1]
        b.getSpring().axis = end.getPos() - start.getPos()
        b.getSpring().pos = start.getPos()  # visualization step
        springStretch = mag(b.getSpring().axis) - b.getConstants()[0]
        force = norm(b.getSpring().axis) * b.getConstants()[1] * springStretch
        # centerofmassVelocity = (start.getVelocity() + end.getVelocity()) / 2
        # start.setVelocity(start.getVelocity() + force * dt / nodeMass - 0.0 * norm(start.getVelocity()) * mag(
        #     start.getVelocity() - centerofmassVelocity) * dt / nodeMass)
        # end.setVelocity(end.getVelocity() - force * dt / nodeMass - 0.0 * norm(end.getVelocity()) * mag(
        #     end.getVelocity() - centerofmassVelocity) * dt / nodeMass)

        centerofmassVelocity = (start.getVelocity() + end.getVelocity()) / 2
        start.setVelocity(start.getVelocity() + force * dt / nodeMass - 0.0 * norm(start.getVelocity()) * mag(
            start.getVelocity() - centerofmassVelocity) * dt / nodeMass)
        end.setVelocity(end.getVelocity() - force * dt / nodeMass - 0.0 * norm(end.getVelocity()) * mag(
            end.getVelocity() - centerofmassVelocity) * dt / nodeMass)

        if springStretch>=0.001:
            b.setConstants([b.getConstants()[0],b.getConstants()[1],b.getConstants()[2],True])
            bonds.remove(b)




    #
    # Image Capture
    #
    # count += 1
    # if count == 1000:
    #     scene.capture("test.png")
    #     count = 0

print("Completed.")
