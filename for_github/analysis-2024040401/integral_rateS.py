import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

# Version 1.0: based on integral_F4.py Version 4.1
#              analyze multiple files - debugged version
#              plot the melting process of multiple initial in one plot.
#              and is used in analysis file
# Add the function to plot the rate of melting at certain instant: melting rate vs time
# Method: Change the analysis_name, analysis_exp, filename_list, C_B, C_list, dt

analysis_name = "analysis-2024040401"
base_path = "/mnt/c/Users/15921/codes/20230831/iceberg-melting-code-master"
path = f"{base_path}/{analysis_name}"
# mode = 0o777
# os.makedirs(path, mode=mode, exist_ok=True)

analysis_exp = "data-456"
filename_list = [
    "data-test2024033122-00", "data-test2024040101-00", "data-test2024040301-00"
]

C_0 = 1
C_B_list = [200, 30, 0]
C_array = np.array(C_B_list) * C_0

dt = 0.5  # the time interval between each data point

imagepath = os.path.join(path, "integral_F.png")

t_matrix = []
t_rate_matrix = []
ice_volume_matrix = []
ice_rate_matrix = []
for filename in tqdm(filename_list):
    filepath = base_path + "/" + filename + "/" + filename + ".h5"

    with h5py.File(filepath, mode='r') as file:
        phase = file['tasks']['f']
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
            matrix_t = matrix_t_array.tolist()
            phase_t.append(matrix_t)

        ice_volume_list = []
        index = 0
        for num in range(len(phase_t)):
            element = phase_t[num]
            # ice_volume_list.append(np.sum(phase_t[num]))
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
    ice_volume_matrix.append(ice_volume_list)

    ice_rate_list = []  # the future melting rate at t
    for index in range(len(ice_volume_list) - 1):
        ice_rate_list.append((-(ice_volume_list[index + 1] - ice_volume_list[index])) / dt)

    ice_rate_matrix.append(ice_rate_list)
    t_matrix.append(t_list)
    t_rate_matrix.append(t_list[:-1])
    # plt.scatter(t_list, ice_volume_list, s=8, label="salinity="+str(C_list[index]))
    # plt.plot(t_list, ice_volume_list, label="salinity="+str(C_list[index]))
    # index += 1

fig, ax = plt.subplots(1, 1, figsize=(10, 6))
# print("length of C_array", len(C_array))
for index in range(len(C_array)):
    ax.plot(t_matrix[index], ice_volume_matrix[index], label=f"salinity = {C_array[index]} g/kg")
ax.set_xlabel('time (s)')
ax.set_ylabel('ice volume (cm^2)')
ax.set_title('ice volume vs time')
ax.legend()
plt.savefig(analysis_name + "_" + analysis_exp + "_F4.png")


fig, ax2 = plt.subplots(1, 1, figsize=(10, 6))
for index in range(len(C_array)):
    ax2.plot(t_rate_matrix[index], ice_rate_matrix[index], label=f"salinity = {C_array[index]} g/kg")
ax2.set_xlabel('time (s)')
ax2.set_ylabel('melting rate (cm^2 / s)')
ax2.set_title('melting rate vs time')
ax2.legend()
plt.savefig(analysis_name + "_" + analysis_exp + "_rate1.png")

plt.show()

plt.close()
