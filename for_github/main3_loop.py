import zescript_3
import merge
import analysis
import seawater as sw

"This code is for the zescript_3 code"
"This is a looped version of main.py that can run several simulations with varying initial conditions"
"Remember to change the FILENAME!!!"

# Change initial conditions here
zescript_3.u_0 = 0


# Change the simulation name and the restart parameter
date = "20240402"
start = 0

zescript_3.restart = start

print("test date is ", date)
flag1 = input("Ready to run the simulation? y/n ")


if flag1 == 'y':
    j = 1
    zescript_3.C_0 = 1
    zescript_3.C_B = 200
    for i in range(2, 22, 2):
        zescript_3.T_B = i
        name = 'test{}{}'.format(date, str(j).zfill(2))
        zescript_3.sim_name = name
        title = 'data-{}-{:0>2d}'.format(name, start)
        merge.sim_name = title
        analysis.filename = title
        zescript_3.main()
        merge.main()
        j += 1

