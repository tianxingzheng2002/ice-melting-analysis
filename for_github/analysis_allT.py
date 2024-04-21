import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
import os

# Updates: converted the units on the salinity scale from C_B percentage to g/kg

"""Important: remember to check the T_B value everytime you run this code !!!!!!"""
"""Important: remember to change the vec_scale_factor, which control the length of the vector field of velocity"""

T_B = 20
vec_scale_factor = 5  # maximum velocity divided by this value gives the velocity represented by the diagonal unit length
vec_interval = 5  # The interval between velocity vectors
vec_color = 'black'
f_min = 0.35
f_max = 0.65  # The range where the water ice interface is identified
f_color = 'white'
f_size = 1

C_B = 200
filename = "data-test2024033122-00"

filepath = filename + "/" + filename + ".h5"
param_name = "parameters-test2024031801-00"
param_path = param_name + "/" + param_name + "_s1/" + param_name + "_s1_p0.h5"

with h5py.File(filepath, mode='r') as file:
    salinity = file['tasks']['T']
    t = salinity.dims[0]['sim_time']
    x = salinity.dims[1][0]
    z = salinity.dims[2][0]
    salinity_list = salinity[:]
    t_list = t[:]
    x_list = x[:]
    z_list = z[:]
    x_v_list = x[::vec_interval]
    z_v_list = z[::vec_interval]
    X_list, Z_list = np.meshgrid(x_list, z_list)
    X_V_list, Z_V_list = np.meshgrid(x_v_list, z_v_list)  # Grids for the velocity vector field.
    salinity_t = []
    for matrix in salinity_list:
        matrix_array = np.array(matrix)
        matrix_t_array = matrix_array.T
        matrix_t = (matrix_t_array * T_B).tolist()
        salinity_t.append(matrix_t[:])

    xvel = file['tasks']['u']
    xvel_list = xvel[:]
    xvel_t = []
    for matrix in xvel_list:
        matrix_array = np.array(matrix)
        matrix_t_array = matrix_array.T
        matrix_t = []
        for line in matrix_t_array:
            matrix_t.append(line[::vec_interval])
        xvel_t.append(matrix_t[::vec_interval])
    xvel_array = np.array(xvel_t)

    zvel = file['tasks']['w']
    zvel_list = zvel[:]
    zvel_t = []
    for matrix in zvel_list:
        matrix_array = np.array(matrix)
        matrix_t_array = matrix_array.T
        matrix_t = []
        for line in matrix_t_array:
            matrix_t.append(line[::vec_interval])
        zvel_t.append(matrix_t[::vec_interval])
    zvel_array = np.array(zvel_t)

    magvel = (xvel_array ** 2 + zvel_array ** 2) ** .5
    lw = (magvel / np.max(magvel))
    magvel_max = np.max(magvel)
    magvel_scale = magvel_max / vec_scale_factor

    phase = file['tasks']['f']
    phase_list = phase[:]
    phase_t = []
    for matrix in phase_list:
        matrix_array = np.array(matrix)
        matrix_t_array = matrix_array.T
        matrix_t = matrix_t_array.tolist()
        phase_t.append(matrix_t)

    ice_volume_list = []
    index = 0
    for num in range(len(phase_t)):
        element = phase_t[num]
        Area = 0
        ice_volume = 0
        for i in range(len(element)):
            for j in range(len(element[i])):
                if j == 0:
                    dx = X_list[i][j + 1] / 2
                elif j == len(element[i]) - 1:
                    dx = (X_list[i][j] - X_list[i][j - 1]) / 2
                else:
                    dx = (X_list[i][j] - X_list[i][j - 1]) / 2 + (X_list[i][j + 1] - X_list[i][j]) / 2
                if i == 0:
                    dz = Z_list[i][j] / 2
                elif i == len(element) - 1:
                    dz = (Z_list[i][j] - Z_list[i - 1][j]) / 2
                else:
                    dz = (Z_list[i][j] - Z_list[i - 1][j]) / 2 + (Z_list[i + 1][j] - Z_list[i][j]) / 2
                Area += dx * dz
                ice_volume += element[i][j] * dx * dz
        ice_volume_list.append(ice_volume)

    x_f_list = []
    z_f_list = []
    phase_plot_list = []
    for t_id in range(len(phase_t)):
        x_f_list.append([])
        z_f_list.append([])
        for i in range(len(phase_t[t_id])):  # i is the index of z_f_list
            for j in range(len(phase_t[t_id][i])):  # j is the index of x_f_list
                if f_max > phase_t[t_id][i][j] > f_min:
                    z_f_list[t_id].append(z_list[i])
                    x_f_list[t_id].append(x_list[j])

    color = "jet"
    cmax = np.max(salinity_t)
    cmin = np.min(salinity_t)
    fig, ax = plt.subplots(figsize=(10, 5))
    image = ax.pcolormesh(X_list, Z_list, salinity_t[0], shading='nearest', cmap=color, vmin=cmin, vmax=cmax)
    scatter = ax.scatter(x_f_list[0], z_f_list[0], color=f_color, s=f_size)
    xvel_half = 0.5 * xvel_array[0]
    zvel_half = 0.5 * zvel_array[0]
    quiver = ax.quiver(X_V_list - xvel_half, Z_V_list - zvel_half, xvel_array[0], zvel_array[0], angles='xy', scale_units='xy', scale=magvel_scale, color=vec_color)
    text = ax.text(0, -0.5, f"time: {round(t_list[0], 1)} s, current ice volume: {round(ice_volume_list[0], 3)} cm^2")
    # print(xvel_array[10])
    print("max velocity is", np.max(magvel))

    def init():
        return [image, scatter, quiver, text]

    def update(frame):
        for coll in ax.collections:
            coll.remove()
        image = ax.pcolormesh(X_list, Z_list, salinity_t[frame], shading='nearest', cmap=color, vmin=cmin, vmax=cmax)
        scatter = ax.scatter(x_f_list[frame], z_f_list[frame], color=f_color, s=f_size)
        u_half = 0.5 * xvel_array[frame]
        w_half = 0.5 * zvel_array[frame]
        quiver = ax.quiver(X_V_list - u_half, Z_V_list - w_half, xvel_array[frame], zvel_array[frame], angles='xy', scale_units='xy', scale=magvel_scale, color=vec_color)
        text.set_text(f"time: {round(t_list[frame], 1)} s, current ice volume: {round(ice_volume_list[frame], 3)} cm^2")
        # image.set_clim(vmin=-2, vmax=2)
        return [image, scatter, quiver, text]

    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title(f'Ice Melting Simulation (velocity field unit diagonal len = {round(magvel_scale, 3)} cm/s) (C0 = {C_B} g/kg)')
    plt.tight_layout()

    animation = FuncAnimation(fig, update, init_func=init, frames=len(salinity_t), interval=100, blit=True)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label('Temperature (degree C)')
    plt.show()

    writervideo = FFMpegWriter(fps=5)
    animation.save(
        "/mnt/c/Users/15921/codes/20230831/iceberg-melting-code-master/" + filename + "/" + filename + "_all_T" + ".mp4",
        writer=writervideo)

    plt.close()
