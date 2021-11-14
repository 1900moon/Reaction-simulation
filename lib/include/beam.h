#pragma once

#include<string>

using namespace std;

class Beam
{
    private:
        double num_cond = 13, num_sigma = 180;
        double input_data[13];
        double cm_angle[180], sigma[180];
        int DIMD = 849;
        double energy1[849],energy2[849];
        double Eloss_d[849],Eloss_C[849];
        double time, intensity, purity, thickness, target_purity, density, strip_x, strip_y, strip_ang, particle_energy;
        double R, detector_sigma;
        double beam_size;
    public:
        void read_parameters();
        void print_cond();
        int get_ini_num();
        void generate_beam(double position[],double direction_Start[], double Inci_Energy[]);
        
        //it is proto type
        void reation_loc_target(double position[], double Inci_Energy[]);
        
        double NuclearReaction(double position[], double direction_Start[], double Inci_Energy[],double Emit_particle[]);
    
        //it is proto type
        double leave_target(double position[],double Emit_particle[],double position_Emit[]);

        int judge_detector(double Emit_particle[], double position_Emit[], double det_position[]);

        double energy_detector(double energy);
};



