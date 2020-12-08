from vpython import *
import numpy as np
import math
# from mainComputation import nodesCoordinates, glassSize
import pandas as pd


#
# Constructing the scene
#
scene2 = canvas(title="Shattering Glass", caption="Animated Display", center=vector(0, 0, 0),
               background=color.black)

floor = box(pos=vector(0, 0, 0), up=norm(vector(0.1, 1, 0)), length=0.4, height=0.001, width=0.1, color=color.green)

t = 0
dt = 0.00001
frame = 30
count = 0
df = pd.read_csv("./output.csv")  # i is the time, and j is the position of nodes
numOfIndices = 27
nodes = []
bondConstants = [0.01 / 3, 99, 0, False]
faces = []
# print(df[("Node " + str(0))][0])


#
# Converting strings to vectors
#
for i in range(numOfIndices):
    temp = df[("Node " + str(i))][0]
    firstComma = temp.find(",")
    secondComma = temp.find(",", firstComma + 1)
    nodes.append(sphere(pos=vector(float(temp[1: firstComma]), float(temp[firstComma + 2: secondComma]), float(temp[secondComma + 2: len(temp) - 1])),
                        radius=bondConstants[0] / 2, color=color.red, visible=False))
    i += 1


while t < 999999999:
    t += 1

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

    for n in range(len(nodes)):
        temp = df[("Node " + str(n))][t]
        firstComma = temp.find(",")
        secondComma = temp.find(",", firstComma + 1)
        nodes[n].pos = vector(float(temp[1: firstComma]), float(temp[firstComma + 2: secondComma]), float(temp[secondComma + 2: len(temp) - 1]))

        #
        # Meshing with triangles
        #
        if t == 20000:
            n1 = 0
            n2 = 0
            maxDistance = 0
            nextDistance = 0
            for j in range(len(nodes)):
                distance = mag(nodes[j].pos - nodes[n].pos)
                if distance > maxDistance:
                    maxDistance = distance
                    n1 = j
                elif distance > nextDistance:
                    nextDistance = distance
                    n2 = j
            T = triangle(
                v0=vertex(pos=nodes[n].pos),
                v1=vertex(pos=nodes[n1].pos),
                v2=vertex(pos=nodes[n2].pos), color=color.red)
            if T in faces:
                T.visible = False
                break
            faces.append(T)

    #
    # Image Capture
    #
    # count += 1
    # if count == trunc(1 / dt / frame):
    #     scene.capture("test.png")
    #     count = 0

print("Visualization Completed.")
