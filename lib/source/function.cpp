#include"../include/function.h"
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

const double STANDARD::HBAR_C = 197.327; //MeV fm
const double STANDARD::ALPHA = 1.0 / 137.036;
const double STANDARD::Qvalue = 7.74939;//MeV
const int STANDARD::DIMD=849;

//unit of mass : MeV/c2
const double MASS::MASS_6HE = 6018885.889 * 1.0e-6 * 931.494013;
const double MASS::MASS_d = 2014101.777844 * 1.0e-6 * 931.494013;
const double MASS::MASS_7Li = 7016003.43426 * 1.0e-6 * 931.494013;
const double MASS::MASS_n = 1008664.91590 * 1.0e-6 * 931.494013;


void test() 
{
    cout << "This is a sample." << endl;
}

double generate_standard() 
{
    random_device rnd; 
    static mt19937 mt(rnd());
    uniform_real_distribution<double> get_rand_uni_real(0.0, 1.0);
    return get_rand_uni_real(mt);
}

double generate_normal(double mu, double sigma) //Box–Muller's method
{
    double value;
    double alpha = generate_standard();
    double beta = generate_standard();
    value = sqrt(-1.0*log(alpha*alpha))*sin(2*M_PI*beta);
    return value*sigma + mu;
}

double Solid_angle(double theta1,double theta2)
{
	double solid_angle;
	double to_rad = M_PI / 180.0;
	theta1 = theta1 * to_rad;
	theta2 = theta2 * to_rad;
	//dΩ=2*Pi*sin(θ)dθ
	solid_angle = 2.0 * M_PI * (cos(theta1)-cos(theta2));
	//cout<<theta1*180.0/M_PI<<"  "<<theta2*180.0/M_PI<<"  "<<solid_angle<<"*****************"<<endl;
	return solid_angle;
}

double Qequation(double Inci_Energy,double theta,double cm_ang_rad, double Gamma)
{
	double mass[4]={MASS::MASS_6HE,MASS::MASS_d,MASS::MASS_7Li,MASS::MASS_n};
	double Qval = STANDARD::Qvalue;
  	double term1,term2,term3,sqrtEb,Eb;
	double theta_cm_turn ;//the emitted partile wiht the largest angle in lab system
 	term1=sqrt(mass[0]*mass[2]*Inci_Energy)*cos(theta)/(mass[3]+mass[2]);
 	term2=((mass[3]-mass[0])/(mass[3]+mass[2])+mass[0]*mass[2]*cos(theta)*cos(theta)/((mass[3]+mass[2])*(mass[3]+mass[2])))*Inci_Energy;
  	term3=mass[3]*Qval/(mass[3]+mass[2]);
 	theta_cm_turn = asin(1.0/Gamma)+90.0*M_PI/180.0;

  	if(cm_ang_rad < theta_cm_turn)
  	{
		sqrtEb = term1+sqrt(abs(term2+term3));
   		Eb = sqrtEb*sqrtEb;
    }
    else
    {
    	sqrtEb = term1-sqrt(abs(term2+term3));
        Eb = sqrtEb*sqrtEb;
    }
  	return Eb;
}
int Dindex(double val)
{
	int ll;
	ll = log(val/10.0)/log(1.02);
	if(ll<0)ll=0;
	if(ll>=STANDARD::DIMD)ll=STANDARD::DIMD-1;
	return ll;
}