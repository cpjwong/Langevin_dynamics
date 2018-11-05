# Langevin_dynamics
Langevin dynamics of a Rouse chain

The program im_bd2.cpp uses implicit Euler method to integrate the Langevin dynamics equation of motion for a Rouse chain (for which the equilibrium bond length is zero.). Rouse mode analysis of the position and velocity of each bead can be obtained using the program, rouse_mode_p2.cpp and rouse_mode_v.cpp, respectively. Visualization of the time correlation functions can be easily achieved using the Python3 programs, such as g_cpp1s, x_corr1.py and plot_sol.py.

The quick and easy way to compile and run all of these programs is by running:

./run.sh
