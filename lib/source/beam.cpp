#include "../include/beam.h"
#include "../include/function.h"
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <cmath>

using namespace std;

void Beam::read_parameters()
{
    ifstream read;
    read.open("/home/zhangq/桌面/C++/Simulation_modified1/condition/input.txt");
    if(!read){
        cout<<"open input.txt error!"<<endl;
        exit(1);
    }
    for(int i=0; i<num_cond; i++){
        read >> input_data[i];
    }
    read.close();
    
    time = input_data[0];
    intensity = input_data[1];
    purity = input_data[2];
    thickness = input_data[3] * 1.0e-4; // cm 
    target_purity = input_data[4];
    density = input_data[5];
    strip_x = input_data[6];
    strip_y = input_data[7];
    strip_ang = input_data[8];
    particle_energy = input_data[9];
    R = input_data[10];
    detector_sigma = input_data[11];
    beam_size = input_data[12];

    //Read differential cross section of the 6He(d,n)7Li reaction
    read.open("/home/zhangq/桌面/C++/Simulation_modified1/condition/sigma.txt");
    if(!read){
        cout<<"Open sigma file error!"<<endl;
        exit(1);
    }
    for(int i =0; i<num_sigma; i++){
        read>>cm_angle[i]>>sigma[i];
    }
    read.close();

	//Read stoping cross section of 7Li in Carbon
    read.open("/home/zhangq/桌面/C++/Simulation_modified1/Stoppingpower/7LiInCarbon.dat");
    if(!read){
        cout<<"Open dEdx1 file error!"<<endl;
        exit(1);
    }
    for(int i =0; i<DIMD; i++){
        read>>energy1[i]>>Eloss_C[i];
    }
    read.close();  

	//Read stoping cross section of 7Li in hydrogen2
    read.open("/home/zhangq/桌面/C++/Simulation_modified1/Stoppingpower/7LiInHydrogen.dat");
    if(!read){
        cout<<"Open dEdx1 file error!"<<endl;
        exit(1);
    }
    for(int i =0; i<DIMD; i++){
        read>>energy2[i]>>Eloss_d[i];
    }
    read.close();  
}

void Beam::print_cond() 
{
    cout << endl;
    cout << "***** initial condition *****" << endl;
    cout << "experiment time  : " << time << " s" << endl;
    cout << "beam intensity : " << intensity << " pps" << endl;
    cout << "purity of main beam : " << purity * 100.0 << " %" << endl;
    cout << "target thickness : " << thickness * 1.0e+4 << " um" << endl;
    cout << "purity of main target : " << target_purity * 100.0 << " %" << endl;
    cout << "target density : " << density << " atoms/cm3" << endl;
    cout << "size of strip : " << strip_x << " x " << strip_y << " cm2" << endl;
    cout << "detector angle : " << strip_ang << " deg" << endl;
    cout << "particle energy : " << particle_energy << " MeV/u" << endl;
    cout << "distance between detector and target : " << R << " cm" << endl;
    cout << "energy resolution : " << detector_sigma << " MeV" << endl;
    cout << "*****************************" << endl;
    cout << endl;
}

// the total number of loop
int Beam::get_ini_num()
{
    return (int)(time  * intensity);
}

void Beam::generate_beam(double position[4],double direction_Start[3],double Inci_Energy[1])
{
    position[0] = beam_size * generate_normal(0.0, 1.0); //cm, size is the radius of the beam spot
    position[1] = beam_size * generate_normal(0.0, 1.0);
    position[2] = -thickness/2.0;
    direction_Start[0] = 0.0;
    direction_Start[1] = 0.0;
    direction_Start[2] = 1.0;
	Inci_Energy[0] = particle_energy;
} 

//Generating the reaction point via random number(***energy loss and straggling from lise++ value ***)
void Beam::reation_loc_target(double position[4], double Inci_Energy[1])
{
    double stop_length = generate_standard() * thickness; //cm
    position[2] += stop_length;
    double energy_loss = 19.891 * stop_length * 1.0e4/1000.0;//energy_loss in unit of MeV, 19.891: keV/um
    double energy_straggling = generate_normal(0.0, 6.31956e-4 * stop_length* 1.0e4 / 50.0);
    Inci_Energy[0] -= energy_loss + energy_straggling;
    position[3] = stop_length;
}

double Beam::NuclearReaction(double position[4],double direction_Start[3], double Inci_Energy[1], double Emit_particle[7])
{
    double mass[4]={MASS::MASS_6HE,MASS::MASS_d,MASS::MASS_7Li,MASS::MASS_n};
	double Qval = STANDARD::Qvalue;
    double delta_cm = 1.0;// 1 degree in the center of mass system
    double solid_angle;
	double Emit_cm_Angle = 180.0;//in the ceter of system
    double to_deg = 180.0 / M_PI;
	double to_rad = M_PI / 180.0;
    double Ecm;
    double cm_ang_deg, cm_ang_rad,cos_lab_angle,lab_angle_rad,lab_angle_deg,omega;
    double cosemitangle = 0.0;
    double Gamma;//cm_to_lab
    double Yield;
    double EmitAngle;
    double EmitEnergy;
	Ecm = Inci_Energy[0] * mass[1]/(mass[0]+mass[1]);

	cm_ang_deg = Emit_cm_Angle * generate_standard();//generate a angle via random number in the c.m. system
    omega=2.0*M_PI*generate_standard();// generate azimuth φ(0-2*PI) via random number
	cm_ang_rad = cm_ang_deg * to_rad;
	
    //calculate the angle in lab system
    Gamma=(mass[0]*mass[2]*(mass[2]+mass[3]))/(mass[1]*mass[3]*(mass[0]+mass[1]))*Ecm/(Ecm+Qval);
    Gamma=sqrt(Gamma);
    cos_lab_angle = (Gamma+cos(cm_ang_rad))/sqrt(1+Gamma*Gamma+2*Gamma*cos(cm_ang_rad));
    lab_angle_rad = acos(cos_lab_angle);
	lab_angle_deg = lab_angle_rad * to_deg;
    //cm_to_lab = pow((1+Gamma*Gamma+2*Gamma*cos(cm_ang_rad)),1.5)/(1+ Gamma*cos(cm_ang_rad));
	
    //calculate the Yield for each angle in the c.m. system
    for(int i=0; i<180; i++){
        if(cm_ang_deg>=i && cm_ang_deg<1){
		    Yield = purity*density*position[3]*target_purity*sigma[i]*1.0e-27*Solid_angle(0.0,cm_angle[i]);
           	//Yield = purity*density*thickness*target_purity*sigma[i]*1.0e-27*Solid_angle(0.0,cm_angle[i]);
            Yield = Yield  * 36.0 * 24.0  * 6.5 * 180;
            break;
        }
        else if( cm_ang_deg == 1.0){
            Yield = purity*density*position[3]*target_purity*sigma[0]*1.0e-27*Solid_angle(cm_angle[0],cm_angle[1]);
            //Yield = purity*density*thickness*target_purity*sigma[0]*1.0e-27*Solid_angle(cm_angle[0],cm_angle[1]);
            Yield = Yield * 36.0 * 24.0  * 6.5 * 180;
            break;
        }
		else if(cm_ang_deg>1 && cm_ang_deg>i && cm_ang_deg<=i+1){
		    Yield = purity*density*position[3]*target_purity*(sigma[i-1]+sigma[i])*0.5*1.0e-27*Solid_angle(cm_angle[i-1],cm_angle[i]);
            //Yield = purity*density*thickness*target_purity*(sigma[i-1]+sigma[i])*0.5*1.0e-27*Solid_angle(cm_angle[i-1],cm_angle[i]);
            Yield = Yield * 36.0 * 24.0  * 6.5 * 180;
            break;
        }
    }

    //Saving the data of directional cosine, angle in lab., Yield, of emitted particle, 
    Emit_particle[0] = sin(lab_angle_rad)*cos(omega);
    Emit_particle[1] = sin(lab_angle_rad)*sin(omega);
    Emit_particle[2] = cos(lab_angle_rad);
    Emit_particle[3] = lab_angle_rad;
    Emit_particle[4] = cm_ang_deg;
    Emit_particle[5] = Yield;

	//calculating the emitted angle and energy of 7Li
	for(int i=0;i<3;i++) cosemitangle=cosemitangle+direction_Start[i]*Emit_particle[i];//the output angle
	EmitAngle = acos(cosemitangle);
	EmitEnergy = Qequation(Inci_Energy[0], EmitAngle, cm_ang_rad, Gamma);
    Emit_particle[6] = EmitEnergy;
	return Yield;
}

//getting the position and energy as 7Li leaves target surface and saving data in a array
//(*** from SRIM value ***) the stopping cross section with ev/(1e15atoms/cm2)
double Beam::leave_target(double position[4], double Emit_particle[7], double position_Emit[3])
{

    double energy_loss;//
    int index;
    index = Dindex(Emit_particle[6]);
    double del_z = thickness/2.0 - position[2];
    double distance = del_z / Emit_particle[2];
    energy_loss = Eloss_d[index]  * 2.0/3.0 + Eloss_C[index]  * 1.0/3.0;//Brag law
    energy_loss = density * energy_loss * 1.0e-6 / (1.0e15) * distance ;//units:MeV/cm
    position_Emit[0] = position[0] + distance * Emit_particle[0];
    position_Emit[1] = position[1] + distance * Emit_particle[1];
    position_Emit[2] = position[2] + distance * Emit_particle[2]; 
    Emit_particle[6] -= energy_loss;
	return 1;
}

//judging whether the 7Li enters in detector surface
int Beam::judge_detector(double Emit_particle[7], double position_Emit[3], double det_position[3])
{
    int flag=0;

    double u = Emit_particle[0];
    double v = Emit_particle[1];
    double w = Emit_particle[2];

    //Rotating the coordinate system around the Y-axis
    double conv_x = position_Emit[0]*cos(strip_ang * M_PI / 180.0) - position_Emit[2]*sin(strip_ang * M_PI / 180.0);
    double conv_y = position_Emit[1];
    double conv_z = position_Emit[2]*cos(strip_ang * M_PI / 180.0) + position_Emit[0]*sin(strip_ang * M_PI / 180.0); 
    
    double conv_u = u*cos(strip_ang * M_PI / 180.0) - w*sin(strip_ang * M_PI / 180.0);
    double conv_v = v;
    double conv_w = w*cos(strip_ang * M_PI / 180.0) + u*sin(strip_ang * M_PI / 180.0);
	
    //flight length as 7Li leaves target
    double factor = (R - conv_z) / conv_w;

    //actual position in detector surface
    conv_x += factor*conv_u;
    conv_y += factor*conv_v;
    
    //judging whether 7Li is in detector
    if(abs(conv_x) < strip_x/2.0 && abs(conv_y) < strip_y/2.0){
        det_position[0] = conv_x;
        det_position[1] = conv_y;
        flag = 1; 
    }
    return flag;
}

double Beam::energy_detector(double energy)
{
    return generate_normal(energy, detector_sigma);
}
