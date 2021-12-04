#!/usr/bin/python3

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#ax.set_xlim3d(-3, 3)
#ax.set_ylim3d(-3, 3)
#ax.set_zlim3d(-3, 3)

# We will be reading files of the form:
#
# 0,0,0
# 0,0,1
# 1,1,0
# 2,1,0
# ...
# (blank)
# HPHPHPHPHP # N/2 letters, where N is the number of given coordinates
#
# If there is no blank space in the end, the plot is made without colors

import csv
from sys import argv

def get_coordinates():
    x = []
    y = []
    z = []
    hp = []

    with open(argv[1]) as fp:
        rd = csv.reader(fp)
        for line in rd:
            if len(line) == 0:
                hp = [i for i in rd.__next__()[0].upper() ]
                break
            x.append(int(line[0])) 
            y.append(int(line[1]))
            z.append(int(line[2]))

    return (x,y,z,hp)


x,y,z,hp = get_coordinates()

hpcolor = 'r'
if len(hp) != 0:
    mapping = {'H': 'r', 'P': 'b'}
    hpcolor = [ mapping[i] for i in hp ]

if len(hp) != 0:
    ax.plot([], [], 'o', c = 'r', label='H')  
    ax.plot([], [], 'o', c = 'b', label='P')  
    ax.legend()

ax.scatter(x[0::2], y[0::2], z[0::2], c='black', marker='x')
ax.scatter(x[1::2], y[1::2], z[1::2], c=hpcolor)

# Plot backbone
ax.plot(x[0::2], y[0::2], z[0::2], c='grey')

# Plot sidechain connections
size = len(x) // 2
for j in range(size):
    i = 2*j
    a = (x[i], x[i+1])
    b = (y[i], y[i+1])
    c = (z[i], z[i+1])
    ax.plot(a, b, c, c='grey')

#for i in range(size):
#    ax.scatter( [i[0] for i in bb], bb[1], bb[2], c='b', marker='^')
#    ax.scatter(sc[0], sc[1], sc[2], c='r')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

if(len(argv) == 3):
    fig.savefig(argv[2])
