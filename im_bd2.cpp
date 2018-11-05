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
//compile with g++ -std=c++11 -Ofast im_bd2.cpp -o im_bd

VectorXd TDMA2(const MatrixXd& A,const VectorXd& b)
{
	int row=A.rows();
	int col=A.cols();
	VectorXd x=VectorXd::Zero(row);
	VectorXd b2=b;
	VectorXd c(row);
	VectorXd a(row);
	for(int i=0;i<row;i++)
	{
		if(i==0)
		{
			c(i)=A(i,i+1);
			a(i)=0;
		}
		else if(i==row-1)
		{
			a(i)=A(i,i-1);
			c(i)=0;
		}
		else
		{
			a(i)=A(i,i-1);
			c(i)=A(i,i+1);
		}
	}
	double err=1;
	double sum1;
	VectorXd d1p(row);
	VectorXd c1p(row);
	VectorXd err1(row);
	int count=0;
	for(int i=0;i<row;i++)
	{
		if(i==0)
		{
			c1p(i)=c(i)/A(i,i);
		}
		else
		{
			c1p(i)=c(i)/(A(i,i)-a(i)*c1p(i-1));
		}
	}
	while(abs(err)>pow(10,-6))
	{
		
		for(int i=0;i<row;i++)
		{
			if(i==0)
			{
				d1p(i)=b2(i)/A(i,i);
			}
			else
			{
				d1p(i)=(b2(i)-a(i)*d1p(i-1))/(A(i,i)-a(i)*c1p(i-1));
			}
		}
		x(row-1)=d1p(row-1);
		for(int i=row-2;i>=0;--i)
		{
			x(i)=d1p(i)-c1p(i)*x(i+1);
		}
		for(int i=0;i<row;i++)
		{
			sum1=0;
			for(int j=0;j<col;j++)
			{
				if(j==i || j==i+1 || j==i-1)
				{
					;
				}
				else
				{
					sum1=sum1+A(i,j)*x(j);
				}
			}
			b2(i)=b(i)-sum1;
		}
		err1=A*x-b;
		err=0;
		for(int i=0;i<row;i++)
		{
			err=err+pow(err1(i),2)/row;
		}
		err=pow(err,0.5);
		count=count+1;
		if(count>1000)
		{
			break;
		}
	}
	return x;
}

int main(int argc, char** argv){
    default_random_engine generator;
	char* filename; 
    int N,t_step,count,n_save,f_save,count2;
    double dt;
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-n") == 0){
			N=stoi(argv[i+1],NULL);//length of polymer
		}
		else if(strcmp(argv[i], "-t") == 0){
			t_step=stoi(argv[i+1],NULL);//t_step
		}
		else if(strcmp(argv[i], "-dt") == 0){
			dt=stod(argv[i+1],NULL);//size of time step dt=0.001 ps
		}
		else if(strcmp(argv[i],"-f")==0){
			filename=argv[i+1];//initial configuration file
		}
		else if(strcmp(argv[i],"-s")==0){
			f_save=stoi(argv[i+1],NULL);//number of steps to be saved
		}
		else{
			;
		}
	}
	n_save=t_step/f_save;
    ifstream infile(filename);
    VectorXd x(N),y(N),z(N),bx(2*N),by(2*N),bz(2*N),rx(2*N),ry(2*N),rz(2*N);
	MatrixXd x1(N,n_save),y1(N,n_save),z1(N,n_save),vx1(N,n_save),vy1(N,n_save),vz1(N,n_save);
	MatrixXd A=MatrixXd::Zero(N,N);
	MatrixXd I=MatrixXd::Identity(N,N);
	MatrixXd I2=MatrixXd::Identity(N,N);
	MatrixXd E=MatrixXd::Identity(N,N);
	MatrixXd Z=MatrixXd::Zero(N,N);
	I2=dt*I2;
    for(int i=0;i<N;i++){
        infile>>rx(i);
        infile>>ry(i);
        infile>>rz(i);
    }
	for(int i=0;i<N;i++){
        rx(i)=0.01*rx(i);
        ry(i)=0.01*ry(i);
        rz(i)=0.01*rz(i);
    }
    double lam,err,rcmx,rcmy,rcmz;
    double k_b=1.38065*pow(10,-27);//A^2 kg ps^-2K^-1
    double T=450;//K
	double xi=5;//ps^-1
    double m=2.327*pow(10,-26);//kg
    double k=1.8*pow(10,-23);//kg ps^-2
	double c2=900;// k/(m)
	double sigma=pow(2*k_b*T*xi*m,0.5);//normal distribution force
	double sigma2=pow(k_b*T/m,0.5);
	E=-(xi*dt)*E;
	normal_distribution<double> distribution(0,sigma);
	normal_distribution<double> distribution2(0,sigma2);
	for(int i=0;i<N;i++){
        rx(i+N)=distribution2(generator);
        ry(i+N)=distribution2(generator);
        rz(i+N)=distribution2(generator);
    }
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	count=0;
	for(int i=0;i<N;i++){
		if(i==0){
			A(i,i)=-c2*dt;
			A(i,i+1)=c2*dt;
		}
		else if(i==N-1){
			A(i,i)=-c2*dt;
			A(i,i-1)=c2*dt;
		}
		else{
			A(i,i)=-2*c2*dt;
			A(i,i-1)=c2*dt;
			A(i,i+1)=c2*dt;
		}
	}
	MatrixXd A1(A.rows()+A.rows(), A.cols()+A.cols());
	A1 << Z, I2,
		  A, E;
	MatrixXd I1=MatrixXd::Identity(2*N,2*N);
	A1=A1-I1;
	for(int it=0;it<t_step;it++){
		if(it==0){
			for(int i=0;i<N;i++){
				x1(i,count)=rx(i);
				y1(i,count)=ry(i);
				z1(i,count)=rz(i);
				vx1(i,count)=rx(i+N);
				vy1(i,count)=ry(i+N);
				vz1(i,count)=rz(i+N);
			}
			count=count+1;
			//cout << "nothing wrong with initialization" << endl;
		}
		else{
			//VectorXd Fr=VectorXd::Zero(3*N);
			for(int i=0;i<N;i++){
				bx(i)=-rx(i);
				by(i)=-ry(i);
				bz(i)=-rz(i);
				bx(i+N)=-rx(i+N)-distribution(generator)*dt/m;
				by(i+N)=-ry(i+N)-distribution(generator)*dt/m;
				bz(i+N)=-rz(i+N)-distribution(generator)*dt/m;
			}
			rx=TDMA2(A1,bx);
			ry=TDMA2(A1,by);
			rz=TDMA2(A1,bz);
			if(remainder(it,f_save)==0){
				for(int i=0;i<N;i++){
					x1(i,count)=rx(i);
					y1(i,count)=ry(i);
					z1(i,count)=rz(i);
					vx1(i,count)=rx(i+N);
					vy1(i,count)=ry(i+N);
					vz1(i,count)=rz(i+N);
				}
				count=count+1;
				cout << "Saved time step= "<< count << endl;
			}
		}
	}
	ofstream myfile;
    myfile.open("trj_x.dat", ios::binary);
    myfile.write((char *) x1.data(), x1.rows() * x1.cols() * sizeof(double));
    myfile.close();
    myfile.open("trj_y.dat", ios::binary);
    myfile.write((char *) y1.data(), y1.rows() * y1.cols() * sizeof(double));
    myfile.close();
    myfile.open("trj_z.dat", ios::binary);
    myfile.write((char *) z1.data(), z1.rows() * z1.cols() * sizeof(double));
    myfile.close();
	myfile.open("trj_vx.dat", ios::binary);
    myfile.write((char *) vx1.data(), vx1.rows() * vx1.cols() * sizeof(double));
    myfile.close();
    myfile.open("trj_vy.dat", ios::binary);
    myfile.write((char *) vy1.data(), vy1.rows() * vy1.cols() * sizeof(double));
    myfile.close();
    myfile.open("trj_vz.dat", ios::binary);
    myfile.write((char *) vz1.data(), vz1.rows() * vz1.cols() * sizeof(double));
    myfile.close();
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    cout << "run time= " << duration << " seconds" << endl;
	return 0;
}



