import zescript_4
import merge
import analysis
import seawater as sw

"This code is for the zescript_4 code"
"This is a looped version of main.py that can run several simulations with varying initial conditions"
"Remember to change the FILENAME!!!"

# Change initial conditions here
zescript_4.u_0 = 0

# Change the simulation name and the restart parameter
date = "20240403"
start = 0

zescript_4.restart = start

print("test date is ", date)
flag1 = input("Ready to run the simulation? y/n ")


if flag1 == 'y':
    j = 2
    zescript_4.C_0 = 1
    zescript_4.T_B = 2
    for i in (200, 30, 0):
        zescript_4.C_B = i
        name = 'test{}{}'.format(date, str(j).zfill(2))
        zescript_4.sim_name = name
        title = 'data-{}-{:0>2d}'.format(name, start)
        merge.sim_name = title
        analysis.filename = title
        zescript_4.main()
        merge.main()
        j += 1
