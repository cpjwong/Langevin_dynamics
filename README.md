# Langevin_dynamics
Langevin dynamics of a Rouse chain

The program im_bd2.cpp uses implicit Euler method to integrate the Langevin dynamics equation of motion for a Rouse chain (for which the equilibrium bond length is zero.). The output files of such program are binary. There should be in total six output files.

Files containing the Cartesian coordinates of the beads: trj_x.dat, trj_y.dat and trj_z.dat
Files containing the velocities of the beads: trj_vx.dat, trj_vy.dat and trj_vz.dat

Rouse mode analysis of the position and velocity of each bead can be obtained using the program, rouse_mode_p2.cpp and rouse_mode_v.cpp, respectively. 

The programs, trj_conv.cpp and trj_conv1.cpp, can be used to chop off the trajectory files of the positions and velocities of the beads, respectively, into any desirable length.

Visualization of the time correlation functions can be easily achieved using the Python3 programs, such as g_cpp1s, x_corr1.py and plot_sol.py.

The quick and easy way to compile and run all of these programs is by running:

./run.sh
