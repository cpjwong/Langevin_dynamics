#include <iostream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <eigen/Eigen/Dense>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <chrono>
#include <string>
#include <omp.h>
using namespace std;
using namespace Eigen;
using namespace std::chrono;
//compile with g++ -std=c++11 -Ofast main2.cpp -o ld4

int main(int argc, char** argv){
	char* filename; 
    int N,t_step,t_step2;
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-n") == 0){
			N=stoi(argv[i+1],NULL);//length of polymer
		}
		else if(strcmp(argv[i], "-ts") == 0){
			t_step=stoi(argv[i+1],NULL);//t_step
		}
		else if(strcmp(argv[i], "-ts2") == 0){
			t_step2=stod(argv[i+1],NULL);//size of time step dt=0.001 ps
		}
		else{
			;
		}
	}
	MatrixXd x(N,t_step);
	MatrixXd y(N,t_step);
	MatrixXd z(N,t_step);
	MatrixXd x2(N,t_step2);
	MatrixXd y2(N,t_step2);
	MatrixXd z2(N,t_step2);
	ifstream infile;
	infile.open("trj_x.dat",ios::binary);	
	infile.read((char *) x.data(),x.rows()*x.cols()*sizeof(double));
	infile.close();
	infile.open("trj_y.dat",ios::binary);
	infile.read((char *) y.data(),y.rows()*y.cols()*sizeof(double));
	infile.close();
	infile.open("trj_z.dat",ios::binary);
	infile.read((char *) z.data(),z.rows()*z.cols()*sizeof(double));
	infile.close();
	int count=0;
	int dt;
	dt=t_step/t_step2;
	cout << dt << endl;
	for(int i=0;i<t_step;i++){
		if(remainder(i,dt)==0){
			for(int i2=0;i2<N;i2++){
				x2(i2,count)=x(i2,i);
				y2(i2,count)=y(i2,i);
				z2(i2,count)=z(i2,i);
				//cout << x(i2,i) << endl;
			}
			count=count+1;
		}
	}
	cout << count << endl;
	ofstream myfile;
    myfile.open("x.dat", ios::binary);
    myfile.write((char *) x2.data(), x2.rows() * x2.cols() * sizeof(double));
    myfile.close();
    myfile.open("y.dat", ios::binary);
    myfile.write((char *) y2.data(), y2.rows() * y2.cols() * sizeof(double));
    myfile.close();
    myfile.open("z.dat", ios::binary);
    myfile.write((char *) z2.data(), z2.rows() * z2.cols() * sizeof(double));
    myfile.close();
	return 0;
}
