import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
import numpy as np
from scipy.spatial import ConvexHull
import sys
import re
 
# Input arguments
# 1 : Input directory
# 2 : Output directory

def parse_points(input, case):
    real = r"\-?\d*\.\d*"
    realex = re.compile(real)
    point = r"\[[(\-?\d*\.\d*), ]*(\-?\d*\.\d*)\]"

    p_r = re.compile(r"[r:]?(" + point + r")")
    m = re.findall(p_r,input)

    direction = [0,0,0]
    xMax = 0.1
    yMax = 0.1
    zMax = 0.1
    pts = []

    i=0
    for i in range(len(m)):
        next = re.findall(realex, m[i][0])
        if i<3:
            for j in range(3):
                if(float(next[j])!=0):
                    direction[j] = float(next[j])
        else:
            x = float(next[0]) - 0.01*direction[0]
            if(x > xMax):
                xMax = x
            y = float(next[1]) - 0.01*direction[1]
            if(y > yMax):
                yMax = y
            z = float(next[2]) - 0.01*direction[2]
            if(z > zMax):
                zMax = 2*z
            pts.append([x,y,z])
        i=i+1

    print("Parsed points for " + case + " :")
    print(pts)
    print("\n")

    for x,y,z in pts:
        if(x < 1):
            xMax = 1
            pts.append([1, y, z])
            if(y>0):
                # xMax = 1
                pts.append([1, 0, z])
                if(z<zMax):
                    pts.append([1, 0, zMax])
            if(z<zMax):
                pts.append([1, y, zMax])
        if(y > 0):
            pts.append([x, 0, z])
            if(z < zMax):
                pts.append([x, 0, zMax])
        if(z < zMax):
            pts.append([x, y, zMax])

    pts = np.array(pts)
    pts = np.unique(pts, axis=0)
    pts[:,0] = -pts[:,0]
    pts[:,2] = -pts[:,2]

    print("Final points for :" + case + " :")
    print(pts)
    print("\n\n")
    return pts, -zMax


def plotConvexHull(ax, pts, case, color):

    # Plotting convex hull
    hull = ConvexHull(pts)
    
    # Plot defining corner points
    ax.plot(pts.T[0], pts.T[1], pts.T[2], color)

    # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
    for s in hull.simplices:
        s = np.append(s, s[0])  # Here we cycle back to the first coordinate
        ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], color)

    # For Showing point values
    # for x, y, z in pts:
    #     print(x, y, z)
    #     label = "(" + str(x) + ", " + str(y) + ", " + str(z) + ")"
    #     ax.text(x, y, z, label)

def plotConvexHull_c_t(ax, pts, case, color):
    hull = ConvexHull(pts)
    ones = -1.2 * np.ones(len(pts))
    ax.plot(ones, pts.T[0], pts.T[1], color)
    for s in hull.simplices:
        on = -1.2 * np.ones(2)
        ax.plot(on, pts[s, 0], pts[s, 1], color)

def plotConvexHull_c_d(ax, pts, case, color):
    hull = ConvexHull(pts)
    ones = 1.2 * np.ones(len(pts))
    ax.plot(pts.T[0], ones, pts.T[1], color)
    for s in hull.simplices:
        on = 1.2 * np.ones(2)
        ax.plot(pts[s, 0], on, pts[s, 1], color)

def plotConvexHull_d_t(ax, pts, case, color):
    hull = ConvexHull(pts)
    ones = zMin * 1.2 * np.ones(len(pts))
    ax.plot(pts.T[0], pts.T[1], ones, color)
    for s in hull.simplices:
        on = zMin * 1.2 * np.ones(2)
        ax.plot(pts[s, 0], pts[s, 1], on, color)


def plot_c_t(ax, pts, case, color):
    new_pts = []
    for p in pts:
        new_pts.append([p[1], p[2]])
    new_pts = np.array(new_pts)
    new_pts = np.unique(new_pts, axis=0)
    plotConvexHull_c_t(ax, new_pts, case, color)

def plot_c_d(ax, pts, case, color):
    new_pts = []
    for p in pts:
        new_pts.append([p[0], p[2]])
    new_pts = np.array(new_pts)
    new_pts = np.unique(new_pts, axis=0)
    plotConvexHull_c_d(ax, new_pts, case, color)

def plot_d_t(ax, pts, case, color):
    new_pts = []
    for p in pts:
        new_pts.append([p[0], p[1]])
    new_pts = np.array(new_pts)
    new_pts = np.unique(new_pts, axis=0)
    plotConvexHull_d_t(ax, new_pts, case, color)

def set_legend(ax):
    ax.plot([10], [10], [100], color = best_col, label = "$G^{+c}$")
    ax.plot([10], [10], [100], color = actual_col, label = "$G$")
    ax.plot([10], [10], [100], color = worst_col, label = "$G^{-c}$")
    # ax.legend()

def get_axes_3d(fig):
    # fig = plt.figure(figsize=(19.8,10.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(-1.2, 0.1)
    ax.set_ylim(-0.1, 1.2)
    ax.set_zlim(zMin*1.2, 0.1)
    ax.set_xlabel('target', fontsize=15)
    ax.set_ylabel('trapped', fontsize=15)
    ax.set_zlabel('cost', fontsize=15)
    set_legend(ax)
    ax.legend()
    ax.set_xticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_zticks([-2, -1.5, -1, -0.5, 0])
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    ax.zaxis.set_tick_params(labelsize=10)
    return ax


def get_axes_c_t(fig):
    # fig = plt.figure(figsize=(19.8,10.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(-1.2, 0)
    ax.set_ylim(0, 1.0)
    ax.set_zlim(zMin, 0.1)
    # ax.set_xlabel('target', fontsize=15)
    ax.set_ylabel('trapped', fontsize=15)
    ax.set_zlabel('cost', fontsize=15, rotation=90)
    set_legend(ax)
    ax.set_xticks([])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_zticks([-2, -1.5, -1, -0.5, 0])
    ax.yaxis.set_tick_params(labelsize=10)
    ax.zaxis.set_tick_params(labelsize=10)
    return ax

def get_axes_c_d(fig):
    # fig = plt.figure(figsize=(19.8,10.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(-1.0, 0)
    ax.set_ylim(0, 1.0)
    ax.set_zlim(zMin, 0)
    ax.set_xlabel('target', fontsize=15)
    # ax.set_ylabel('trapped', fontsize=15)
    ax.set_zlabel('cost', fontsize=15, rotation=90)
    set_legend(ax)
    ax.set_xticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0])
    ax.set_yticks([])
    ax.set_zticks([-2, -1.5, -1, -0.5, 0])
    ax.xaxis.set_tick_params(labelsize=10)
    ax.zaxis.set_tick_params(labelsize=10)
    return ax

def get_axes_d_t(fig):
    # fig = plt.figure(figsize=(19.8,10.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(-1.0, 0)
    ax.set_ylim(-0, 1.0)
    ax.set_zlim(zMin, 0)
    ax.set_xlabel('target', fontsize=15)
    ax.set_ylabel('trapped', fontsize=15)
    # ax.set_zlabel('cost', fontsize=15)
    set_legend(ax)
    ax.set_xticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_zticks([])
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    return ax


model_f = open(sys.argv[1] + "model_pareto.txt")
best_f = open(sys.argv[1] + "best_pareto.txt")
worst_f = open(sys.argv[1] + "worst_pareto.txt")

model_input = model_f.read()
best_input = best_f.read()
worst_input = worst_f.read()

zMin = 0.1
pts_model, zMin1 = parse_points(model_input, "model")
pts_best, zMin2 = parse_points(best_input, "best")
pts_worst, zMin3 = parse_points(worst_input, "worst")
zMin = min(zMin, zMin1, zMin2, zMin3)

fig = plt.figure(figsize=(19.8,10.0))
# fig.SubplotParams(left=0, right=0.5, bottom=0, top=0.5)

actual_col = "blue"
best_col = "#6eb84f"
worst_col = "#db6244"


ax1 = get_axes_c_d(fig)
# tmp_planes = ax1.zaxis._PLANES 
# ax.zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], 
#                      tmp_planes[0], tmp_planes[1], 
#                      tmp_planes[4], tmp_planes[5])
plot_c_d(ax1, pts_model, "model", actual_col)
plot_c_d(ax1, pts_best, "best", best_col)
plot_c_d(ax1, pts_worst, "worst", worst_col)
ax1.view_init(elev=0.01, azim=270.01)
plt.savefig(sys.argv[2] + "c_vs_d.png", format = 'png', dpi = 300, bbox_inches='tight')

ax2 = get_axes_c_t(fig)
plot_c_t(ax2, pts_model, "model", actual_col)
plot_c_t(ax2, pts_best, "best", best_col)
plot_c_t(ax2, pts_worst, "worst", worst_col)
ax2.view_init(elev=0.01, azim=-0.01)
plt.savefig(sys.argv[2] + "c_vs_t.png", format = 'png', dpi = 300, bbox_inches='tight')

ax3 = get_axes_d_t(fig)
plot_d_t(ax3, pts_model, "model", actual_col)
plot_d_t(ax3, pts_best, "best", best_col)
plot_d_t(ax3, pts_worst, "worst", worst_col)
ax3.view_init(elev=89.99, azim=270)
plt.savefig(sys.argv[2] + "d_vs_t.png", format = 'png', dpi = 300, bbox_inches='tight')

ax = get_axes_3d(fig)
plotConvexHull(ax, pts_model, "model", actual_col)
plotConvexHull(ax, pts_best, "best", best_col)
plotConvexHull(ax, pts_worst, "worst", worst_col)
ax.view_init(elev=15, azim=-45)
plt.savefig(sys.argv[2] + "3d.png", format = 'png', dpi = 300, bbox_inches='tight')

# handles, labels = ax.get_legend_handles_labels()
# fig.legend(handles, labels, loc='center')
# plt.savefig(sys.argv[2] + "plot.png", format = 'png', dpi = 300, bbox_inches='tight')
# plt.show()

