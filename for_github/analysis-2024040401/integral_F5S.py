import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

# Version 5.0: analyze multiple files - in individual analysis file folders and only plot end state
#              plot end-state volume vs initial salinity
# Version 5.0 Super: analyze rows of experiments with different initial temperatures
# Method: Change the analysis_name, filename_list, C_B, C_list

analysis_name = "analysis-2024040401"
base_path = "/mnt/c/Users/15921/codes/20230831/iceberg-melting-code-master"
path = f"{base_path}/{analysis_name}"


analysis_exp = "data-ABCD"
T_B1 = 2
filename_list1 = [
    "data-test2024033001-00", "data-test2024033002-00", "data-test2024033003-00",
    "data-test2024033004-00", "data-test2024033005-00", "data-test2024033006-00",
    "data-test2024033007-00", "data-test2024033112-00", "data-test2024033008-00", "data-test2024033009-00",
    "data-test2024033010-00", "data-test2024033011-00", "data-test2024033113-00",
    "data-test2024033114-00", "data-test2024033115-00", "data-test2024033116-00",
    "data-test2024033117-00", "data-test2024033118-00", "data-test2024033119-00",
    "data-test2024033120-00", "data-test2024033121-00"
]

T_B2 = 20
filename_list2 = [
    "data-test2024033101-00", "data-test2024033102-00", "data-test2024033103-00",
    "data-test2024033104-00", "data-test2024033105-00", "data-test2024033106-00",
    "data-test2024033107-00", "data-test2024032402-00", "data-test2024033108-00", "data-test2024033109-00",
    "data-test2024033110-00", "data-test2024033111-00",
    "data-test2024032403-00",
    "data-test2024032404-00", "data-test2024032405-00", "data-test2024032406-00",
    "data-test2024032407-00", "data-test2024032408-00", "data-test2024032409-00",
    "data-test2024032410-00", "data-test2024032411-00"
]

filename_list = filename_list1 + filename_list2

C_0 = 1
C_B_list = [0, 3, 6, 9, 12, 15, 18, 20, 21, 24, 27, 30, 40, 60, 80, 100, 120, 140, 160, 180, 200]
C_array = np.array(C_B_list) * C_0

fig, ax = plt.subplots(figsize=(10, 6.5))
ax.set_title('Remaining Ice Volume after 25s vs Initial Salinity')

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

plt.scatter(C_array, ice_volume_end_list[:21], s=8, c='b', label=f"initial temperature = {T_B1} degree C")
plt.plot(C_array, ice_volume_end_list[:21], c='b')
plt.scatter(C_array, ice_volume_end_list[21:], s=8, c='red', label=f"initial temperature = {T_B2} degree C")
plt.plot(C_array, ice_volume_end_list[21:], c='red')
plt.xlabel('Salinity (g/kg)')
plt.ylabel('ice volume (cm^2)')
plt.legend()
plt.savefig(analysis_name + "_" + analysis_exp + "_F5.png")
plt.show()
