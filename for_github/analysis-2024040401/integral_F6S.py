import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

# Version 6.0: analyze multiple files - in individual analysis file folders and only plot end state
#              plot end-state ice volume vs initial temperature
# Version 6.0 Super: analyze multiple salinity
# Method: Change the analysis_name, filename_list, T_B, analysis_exp, C_B

analysis_name = "analysis-2024040401"
base_path = "/mnt/c/Users/15921/codes/20230831/iceberg-melting-code-master"
path = f"{base_path}/{analysis_name}"

# Experiment E: Uncomment the following to analyze E column
analysis_exp1 = "data-E"
filename_list1 = [
    "data-test2024032501-00", "data-test2024032502-00", "data-test2024032503-00",
    "data-test2024032504-00", "data-test2024032505-00", "data-test2024032506-00",
    "data-test2024032507-00", "data-test2024032508-00", "data-test2024032509-00",
    "data-test2024032510-00"
]
C_B1 = 0

# Experiment F: Uncomment the following to analyze F column
analysis_exp2 = "data-F"
filename_list2 = [
    "data-test2024040102-00", "data-test2024040103-00",
    "data-test2024040104-00", "data-test2024040105-00", "data-test2024040106-00",
    "data-test2024040107-00", "data-test2024040108-00", "data-test2024040109-00",
    "data-test2024040110-00", "data-test2024040111-00"
]
C_B2 = 30

# Experiment G: Uncomment the following to analyze G column
analysis_exp3 = "data-G"
filename_list3 = [
    "data-test2024040201-00", "data-test2024040202-00", "data-test2024040203-00",
    "data-test2024040204-00", "data-test2024040205-00", "data-test2024040206-00",
    "data-test2024040207-00", "data-test2024040208-00", "data-test2024040209-00",
    "data-test2024040210-00"
]
C_B3 = 200

T_B_list = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
T_array = np.array(T_B_list)
filename_list = filename_list1 + filename_list2 + filename_list3
analysis_exp = "data-EFG"


fig, ax = plt.subplots(figsize=(10, 6.5))
ax.set_title('Remaining Ice Volume after 25s vs Initial Temperature')

imagepath = os.path.join(path, "integral_F.png")

t_matrix = []
ice_volume_matrix = []
ice_volume_end_list = []
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

        element = phase_t[-1]
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
        ice_volume_end_list.append(ice_volume)
    t_matrix.append(t_list)
    # plt.scatter(t_list, ice_volume_list, s=8, label="salinity="+str(C_list[index]))
    # plt.plot(t_list, ice_volume_list, label="salinity="+str(C_list[index]))
    # index += 1

plt.scatter(T_array, ice_volume_end_list[:10], s=8, c='aqua', label=f"initial salinity = {C_B1} g/kg")
plt.plot(T_array, ice_volume_end_list[:10], c='aqua')
plt.scatter(T_array, ice_volume_end_list[10:20], s=8, c='royalblue', label=f"initial salinity = {C_B2} g/kg")
plt.plot(T_array, ice_volume_end_list[10:20], c='royalblue')
plt.scatter(T_array, ice_volume_end_list[20:], s=8, c='purple', label=f"initial salinity = {C_B3} g/kg")
plt.plot(T_array, ice_volume_end_list[20:], c='purple')
plt.xlabel('Temperature (degree C)')
plt.ylabel('ice volume (cm^2)')
plt.legend()
plt.savefig(analysis_name + "_" + analysis_exp + "_F6.png")
plt.show()
