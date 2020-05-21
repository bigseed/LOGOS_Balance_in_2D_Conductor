import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)


def make_vertex():
    vertex = list()
    for theta in np.arange(0, 2*np.pi, 2*np.pi/n):
        new_x, new_y = rotate(scale*np.cos(np.pi/n),
                              scale*np.sin(np.pi/n), theta)
        vertex.append((new_x, new_y))
    return vertex


def rotate(x, y, theta):
    rotate_x = x*np.cos(theta) - y*np.sin(theta)
    rotate_y = x*np.sin(theta) + y*np.cos(theta)

    return rotate_x, rotate_y


def spread_in_part():
    x = np.random.uniform(0, scale*np.cos(np.pi/n), number_of_pts)
    y = np.zeros(number_of_pts)

    for idx, val in enumerate(x):
        y[idx] = np.random.uniform(-val *
                                   np.tan(np.pi/n), val*np.tan(np.pi/n))

    x, y = np.reshape(x, (number_of_pts, 1)), np.reshape(y, (number_of_pts, 1))
    return (x, y)


def spread_pts():
    x, y = spread_in_part()
    electrons = np.hstack([x, y])
    for i in range(1, n):
        x, y = spread_in_part()
        x, y = rotate(x, y, 2*np.pi/n*i)
        temp_electrons = np.hstack([x, y])
        electrons = np.vstack([electrons, temp_electrons])

    return electrons


def is_out(x, y):
    angle_subtracted = 0

    while x / distance(x, y, 0, 0) < np.cos(np.pi/n):
        x, y = rotate(x, y, -2*np.pi/n)
        angle_subtracted += 2*np.pi/n

    if x > scale * np.cos(np.pi/n):
        tan = y / x
        x = scale * np.cos(np.pi/n)
        y = tan * x

        x, y = rotate(x, y, angle_subtracted)

        return (True, (x, y))
    else:
        return (False, 0)


def calc_force(e1, e2):
    k = 8.987e9
    x1, y1, x2, y2 = e1[0], e1[1], e2[0], e2[1]

    if distance(x1, y1, x2, y2) != 0:
        F = k / distance(x1, y1, x2, y2)**3
        fx1, fy1 = F*(x1 - x2), F*(y1 - y2)
        fx2, fy2 = F*(x2 - x1), F*(y2 - y1)
        return (fx1, fy1, fx2, fy2)
    return (0, 0, 0, 0)


n = 3
scale = 15
number_of_pts = 150
me = 9.109e-31
iteration = 10
heatamp_scale = 1/4

vertex = make_vertex()
pts = spread_pts()
pts = np.hstack([pts, np.zeros((n*number_of_pts, 2))])

for i in range(iteration):
    print(f'Iter {i + 1}')

    for idx_1, e_1 in enumerate(pts):
        for idx_2, e_2 in enumerate(pts):
            if idx_2 <= idx_1:
                continue
            fx1, fy1, fx2, fy2 = calc_force(e_1, e_2)
            e_1[2] += fx1
            e_1[3] += fy1
            e_2[2] += fx2
            e_2[3] += fy2

    for electron in pts:
        electron[0] += electron[2]
        electron[1] += electron[3]

        is_outside = is_out(electron[0], electron[1])

        if is_outside[0]:
            electron[0], electron[1] = is_outside[1][0], is_outside[1][1]

    pts[:, 2:] = 0


if pd.isnull(pts).any():
    print(f'There is {pd.isnull(pts).sum()} null!')

if True:
    pts = np.around(pts/heatamp_scale).astype('int')
    grid_length = int(scale/heatamp_scale)
    how_many_pts = np.zeros((2*grid_length + 1, 2*grid_length + 1))
    for electron in pts:
        x, y = electron[0], electron[1]
        how_many_pts[y + grid_length][x + grid_length] += 1

    df = pd.DataFrame(how_many_pts, index=np.arange(-grid_length, grid_length + 1),
                      columns=np.arange(-grid_length, grid_length + 1))

    ax = sns.heatmap(df)

    plt.title(f'n{n}_iter{iteration}')
    plt.savefig(f'n{n}_iter{iteration}.png')

if not True:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for i in range(n):
        ax.plot((vertex[i-1][0], vertex[i][0]),
                (vertex[i-1][1], vertex[i][1]), c='red')
    ax.set_xlim(-scale*1.2, scale*1.2)
    ax.set_ylim(-scale*1.2, scale*1.2)
    ax.scatter(pts[:, 0], pts[:, 1], s=0.5)
    plt.title(f'n{n}_electron_randomly_spreaded')
    plt.savefig(f'n{n}.png')
