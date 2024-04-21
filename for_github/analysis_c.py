import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
import os

# Updates: converted the units on the salinity scale from C_B percentage to g/kg

"""Important: remember to check the C_B value everytime you run this code !!!!!!"""

C_B = 20

filename = "data-test2024022511-00"
filepath = filename + "/" + filename + ".h5"
param_name = "parameters-test2024012101-00"
param_path = param_name + "/" + param_name + "_s1/" + param_name + "_s1_p0.h5"

# with h5py.File(param_path, mode='r') as file:
#     Sc = file['tasks']['Sc']
#     Sc_num = Sc[:][0][0][0]
#
# print(Sc_num)

with h5py.File(filepath, mode='r') as file:
    phase = file['tasks']['C']
    t = phase.dims[0]['sim_time']
    x = phase.dims[1][0]
    z = phase.dims[2][0]
    phase_list = phase[:]
    t_list = t[:]
    x_list = x[:]
    z_list = z[:]
    X_list, Z_list = np.meshgrid(x_list, z_list)
    phase_t = []
    for matrix in phase_list:
        matrix_array = np.array(matrix)
        matrix_t_array = matrix_array.T
        matrix_t = (matrix_t_array * C_B).tolist()
        phase_t.append(matrix_t)
    test = np.random.rand(256, 128)
    color = "jet"
    cmax = np.max(phase_t)
    cmin = np.min(phase_t)

    # cmax = 0.1
    # cmin = 0.04
    fig, ax = plt.subplots(figsize=(10, 5))
    image = ax.pcolormesh(X_list, Z_list, phase_t[0], shading='nearest', cmap=color, vmin=cmin, vmax=cmax)
    # image.set_clim(vmin=-2, vmax=2)

    def init():
        return [image]

    def update(frame):
        for coll in ax.collections:
            coll.remove()

        image = ax.pcolormesh(X_list, Z_list, phase_t[frame], shading='nearest', cmap=color, vmin=cmin, vmax=cmax)
        # image.set_clim(vmin=-2, vmax=2)
        return [image]

    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title('Ice Melting Simulation')
    plt.tight_layout()

    animation = FuncAnimation(fig, update, init_func=init, frames=len(phase_t), interval=50, blit=True)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label('Salinity (g/kg)')
    plt.show()

    writervideo = FFMpegWriter(fps=10)
    animation.save(
        "/mnt/c/Users/15921/codes/20230831/iceberg-melting-code-master/" + filename + "/" + filename + "_C" + ".mp4",
        writer=writervideo)

    plt.close()