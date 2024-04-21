import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import seawater as sw

# Version 1.0: plot the total kinetic energy change: total KE vs time
# Method: Change the analysis_name, analysis_exp, filename_list, C_B_list, C_0,
#         rho0_list (the codes for calculating rho0), T_B

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

T_B = 20

rho0_list = []
for C_real in C_array:
    rho0_list.append(sw.dens0(s=C_real, t=T_B))

fig, ax = plt.subplots(figsize=(10, 6))
ax.set_title('Total Kinetic Energy vs Time')

imagepath = os.path.join(path, "integral_F.png")

t_matrix = []
ke_matrix = []
file_id = 0
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

        xvel = file['tasks']['u']
        xvel_list = xvel[:]
        xvel_t = []
        for matrix in xvel_list:
            matrix_array = np.array(matrix)
            matrix_t_array = matrix_array.T
            matrix_t = matrix_t_array.tolist()
            xvel_t.append(matrix_t)

        zvel = file['tasks']['w']
        zvel_list = zvel[:]
        zvel_t = []
        for matrix in zvel_list:
            matrix_array = np.array(matrix)
            matrix_t_array = matrix_array.T
            matrix_t = matrix_t_array.tolist()
            zvel_t.append(matrix_t)

        ke_list = []  # kinetic energy list
        index = 0
        for num in range(len(phase_t)):
            element = phase_t[num]
            u_matrix = xvel_t[num]
            w_matrix = zvel_t[num]
            # ice_volume_list.append(np.sum(phase_t[num]))
            Area = 0
            ke_total = 0
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
                    ke_total += 0.5 * (u_matrix[i][j]**2 + w_matrix[i][j]**2) * rho0_list[file_id] * dx * dz
            ke_list.append(ke_total)
    ke_matrix.append(ke_list)
    t_matrix.append(t_list)
    file_id += 1
    # plt.scatter(t_list, ice_volume_list, s=8, label="salinity="+str(C_list[index]))
    # plt.plot(t_list, ice_volume_list, label="salinity="+str(C_list[index]))
    # index += 1

label_list = ["rho0 at s=200g/kg", "rho0 at s=5g/kg"]
for index in range(len(C_array)):
    plt.plot(t_matrix[index], ke_matrix[index], label=f"salinity = {C_array[index]} g/kg")
plt.xlabel('time (s)')
plt.ylabel('total kinetic energy (10^-10 J)')
plt.legend()
plt.savefig(analysis_name + "_" + analysis_exp + "_KE1.png")
plt.show()
