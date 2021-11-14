#pragma once

#include<string>

using namespace std;

struct STANDARD
{
  static const double HBAR_C;
  static const double ALPHA;
  static const double Qvalue;
  static const int DIMD;
};

struct MASS
{
  static const double MASS_6HE;
  static const double MASS_3H;

  static const double MASS_d;
  static const double MASS_12C;

  static const double MASS_7Li;
  static const double MASS_n;
};

void test();

double generate_standard();
double generate_normal(double mu, double sigma);
double Solid_angle(double theta1,double theta2);
double Qequation(double Inci_Energy,double theta,double cm_ang_rad, double Gamma);
int Dindex(double val);