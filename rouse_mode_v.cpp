#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <map>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <eigen/Eigen/Dense>
#include <chrono>
using namespace Eigen;
using namespace std;
using namespace std::chrono;
double A_f(int p,int n,int N){
	return cos((double)(p)*n*M_PI/(double)(N));
}

int main(int argc, char** argv)
{
	cout << "For rouse analysis with small time step." << endl;
	default_random_engine generator;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int N,t_step,an,cn,ac;
	char* filename;
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i], "-ts") == 0){
			t_step=stoi(argv[i+1],NULL);
		}
		else if(strcmp(argv[i],"-an")==0){
			an=stoi(argv[i+1],NULL);
		}
		else if(strcmp(argv[i],"-cn")==0){
			cn=stoi(argv[i+1],NULL);
		}
		else{
			;
		}
	}
	ac=an/cn;
	MatrixXd x(an,t_step);
	MatrixXd y(an,t_step);
	MatrixXd z(an,t_step);
	MatrixXd vx(an,t_step);
	MatrixXd vy(an,t_step);
	MatrixXd vz(an,t_step);
	ifstream infile;
	infile.open("vx.dat",ios::binary);
	infile.read((char *) vx.data(),vx.rows()*vx.cols()*sizeof(double));
	infile.close();
	infile.open("vy.dat",ios::binary);
	infile.read((char *) vy.data(),vy.rows()*vy.cols()*sizeof(double));
	infile.close();
	infile.open("vz.dat",ios::binary);
	infile.read((char *) vz.data(),vz.rows()*vz.cols()*sizeof(double));
	infile.close();
	cout << "Reading data is complete!" << endl;
	cout << "first element of vx= " << vx(0,0) << endl;
	cout << "last element of vx= " << vx(an-1,t_step-1) << endl;
	cout << ac << endl;
	MatrixXd Xp1(ac,t_step),Xp2(ac,t_step),Xp3(ac,t_step);
	int tau_size;
	double c;
	double sum1,sum2,sum3;
	for(int i=0;i<cn;i++){
		for(int j=0;j<t_step;j++){
			for(int k1=0;k1<ac;k1++){
				sum1=0;
				sum2=0;
				sum3=0;
				for(int k2=0;k2<ac;k2++){
					sum1=sum1+vx(k2+i*ac,j)*A_f(k1,k2,ac);
					sum2=sum2+vy(k2+i*ac,j)*A_f(k1,k2,ac);
					sum3=sum3+vz(k2+i*ac,j)*A_f(k1,k2,ac);
				}
				Xp1(k1,j)=sum1/(double)ac;
				Xp2(k1,j)=sum2/(double)ac;
				Xp3(k1,j)=sum3/(double)ac;
			}
		}
		cout << "Conversion to normal coordinate is complete! " << endl;
		VectorXd mean_Xsq=VectorXd::Zero(ac);
		for(int j=0;j<ac;j++){
			for(int k=0;k<t_step;k++){
				mean_Xsq(j)=mean_Xsq(j)+(Xp1(j,k)*Xp1(j,k)+Xp2(j,k)*Xp2(j,k)+Xp3(j,k)*Xp3(j,k))/(double)t_step;
			}			
		}
		MatrixXd corr(t_step,ac);
		for(int k=0;k<ac;k++){
			for(int j=0;j<t_step;j++){
				if(j == 0){
					corr(j,k)=1;
				}
				else{
					tau_size=t_step-j;
					c=0;
					for (int k1=0;k1<tau_size;k1++){
						c=c+(Xp1(k,k1)*Xp1(k,k1+j)+Xp2(k,k1)*Xp2(k,k1+j)+Xp3(k,k1)*Xp3(k,k1+j))/(mean_Xsq(k)*(double)(tau_size));
					}
					corr(j,k)=c;
				}
			}
		}
		cout << "Correlation function calculation complete! chain number= " << i << endl;
		string file1="Vp1_";
		string file2="Vp2_";
		string file3="Vp3_";
		string file4="rousev_1_";
		string file5="mean_sq_V_";
		ofstream myfile;
		myfile.open((file1+to_string(i)).c_str(), ios::binary);
		myfile.write((char *) Xp1.data(), Xp1.rows() * Xp1.cols() * sizeof(double));
		myfile.close();
		myfile.open((file2+to_string(i)).c_str(), ios::binary);
		myfile.write((char *) Xp2.data(), Xp2.rows() * Xp2.cols() * sizeof(double));
		myfile.close();
		myfile.open((file3+to_string(i)).c_str(), ios::binary);
		myfile.write((char *) Xp3.data(), Xp3.rows() * Xp3.cols() * sizeof(double));
		myfile.close();
		myfile.open((file4+to_string(i)).c_str(), ios::binary);
		myfile.write((char *) corr.data(), corr.rows() * corr.cols() * sizeof(double));
		myfile.close();
		myfile.open((file5+to_string(i)).c_str(), ios::binary);
		myfile.write((char *) mean_Xsq.data(), mean_Xsq.rows() * sizeof(double));
		myfile.close();
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    cout << "run time= " << duration << " seconds" << endl;
	return 0;
}
