#!/bin/bash

g++ -std=c++11 -Ofast -fopenmp im_bd2.cpp -o im_ld
g++ -std=c++11 -Ofast rouse_mode_p2.cpp -o rmode_ld_linear
g++ -std=c++11 -Ofast rouse_mode_v.cpp -o rmode_v
g++ -std=c++11 -Ofast trj_conv.cpp -o trj_conv
g++ -std=c++11 -Ofast trj_conv1.cpp -o trj_conv_vel
echo "Compilation complete!"
time=3500
ts1=0.01
./im_ld -f initial2.txt -n 50 -t $time -dt $ts1 -s 1
./trj_conv_vel -n 50 -ts $time -ts2 $time
./rmode_v -ts $time -cn 1 -an 50
python3 g_cpp1v.py -cn 1 -an 50 -t $ts1 -s $time
./trj_conv -n 50 -ts $time -ts2 $time
./rmode_ld_linear -ts $time -cn 1 -an 50
python3 g_cpp1s -cn 1 -an 50 -t $ts1 -s $time
