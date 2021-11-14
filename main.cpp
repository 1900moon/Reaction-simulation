#include <function.h>
#include <beam.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main( int argc, char **argv )
{

    //prepare for rootfile
    double ini_x, ini_y, ini_z, ini_energy;//Initial incident particle information
    double reac_x, reac_y, reac_z, reac_energy;// The reaction point information
    double Yield_reac,Yield;//The reaction yield of the 6He(d,n)7Li reaction
    double cm_ang;//The random emit angle in the c.m. system
    double part_energy, part_det_energy;//the energy of 7Li left the target surface, and the detected 7Li with Si detector
    double det_x, det_y;// the 7Li position in Detector surface
   
    //Save project information in TREE
    TTree *tree = new TTree("tree", "elastic scattering");

    //Initial particle information
    tree->Branch("ini_x", &ini_x, "ini_x/D");
    tree->Branch("ini_y", &ini_y, "ini_y/D");
    tree->Branch("ini_z", &ini_z, "ini_z/D");
    tree->Branch("ini_energy", &ini_energy, "ini_energy/D");

    //The reaction point information:position, reaction energy, 
    tree->Branch("reac_x", &reac_x, "reac_x/D");
    tree->Branch("reac_y", &reac_y, "reac_y/D");
    tree->Branch("reac_z", &reac_z, "reac_z/D");
    tree->Branch("reac_energy", &reac_energy, "reac_energy/D");
    
	//The emmited particle of 7Li:emited angle in the c.m. system , Yield of the d(6He,7Li)n
    tree->Branch("cm_ang", &cm_ang, "cm_ang/D");
    tree->Branch("Yield", &Yield, "Yield/D");

	//The emitted energy from target surface and the detected energy within the si detector
    tree->Branch("part_energy", &part_energy, "part_energy/D");
    tree->Branch("part_det_energy", &part_det_energy, "part_det_energy/D");
    //The 7Li particle position in detector surface
   
    tree->Branch("det_x", &det_x, "det_x/D");
    tree->Branch("det_y", &det_y, "det_y/D");


    //Set primary beam and read parameters
    Beam *beam_test = new Beam();
    beam_test->read_parameters();
    beam_test->print_cond();

    //position[0]: location_x (cm)
    //position[1]: location_y (cm)
    //position[2]: location_z (cm)
    //position[3: Flightlength
    //Inci_energy[0]:the initial energy
    double position[4],direction_Start[3];//the incident particle's position and direction
    double Emit_particle[7],position_Emit[3];// the emitted 7Li's direction and position
    double det_position[3];//the emitted 7Li position in detector surface
    double Inci_Energy[1];//the initial energy of the incident particle 6He


    double Yield_all[180];
    int cm_all[180];
	for(int i=0; i<180; i++){
		cm_all[i] = 0.0;
		Yield_all[i] = 0.0;
    }

    TH1F *h_strip = new TH1F( "h_strip", "", 100, 0., 10 );

    //Start simulation
    for(int loop=0; loop<beam_test->get_ini_num(); loop++){
        ini_x=0.0, ini_y=0.0, ini_z=0.0, ini_energy=0.0;
        reac_x=0.0, reac_y=0.0, reac_z=0.0, reac_energy=0.0;
        cm_ang=0.0, Yield=0.0;
        part_energy=0.0, part_det_energy=0.0;
        det_x=0.0, det_y=0.0;    
        if((loop+1)%10==0){
            cout << "\r" << "proceeding... " << 100.0 * (double)(loop+1)/(double)beam_test->get_ini_num() << " %" << string(20, ' ');
        }

        // get the initial information of incident particle
        beam_test->generate_beam(position,direction_Start,Inci_Energy);      
        ini_x = position[0];
        ini_y = position[1];
        ini_z = position[2];
        ini_energy = Inci_Energy[0];

        //get the reaction point information including position and energy slowed by target
        beam_test->reation_loc_target(position,Inci_Energy);
		reac_x = position[0];
		reac_y = position[1];
		reac_z = position[2];
		reac_energy = Inci_Energy[0];
		
        //get the reaction Yield and direction of emitted particle 7Li
		Yield_reac = beam_test->NuclearReaction(position,direction_Start, Inci_Energy, Emit_particle);

        //get the position and energy as 7Li leaves the target surface
		beam_test->leave_target(position, Emit_particle, position_Emit);
		part_energy = Emit_particle[6]/7.0;
		
        // judging whether 7Li is in detector
        int flag;
        flag = beam_test->judge_detector(Emit_particle, position_Emit, det_position);
		
        //collecting the information of 7Li when they enter in detector surface
        if(flag == 1){
            det_x = det_position[0];
            det_y = det_position[1];
            cm_ang = Emit_particle[4];
            Yield = Yield_reac;
			
            //count the total Yield for each angle(1-180 degrees)
            for(int i=0;i<180;i++)
            {
                if(cm_ang>=i && cm_ang< i+1)
                {
                    cm_all[i]=i+1;
                    Yield_all[i]= Yield_all[i] + Yield;
                    break;
                }
                if(cm_ang == 180.0)
                {
                    Yield_all[179]=Yield_all[179] + Yield;
                    break;
                }
            }
            part_det_energy = beam_test->energy_detector(Emit_particle[6]);
            part_det_energy = part_det_energy/7.0;
            h_strip->Fill(part_det_energy);
            tree->Fill();
        }
       // tree->Fill();
    }
    ofstream read;
    read.open("Yield.txt");
    for(int i=0;i<180;i++){
        cout<<cm_all[i]<<"  "<<Yield_all[i]<<"  "<<endl;
        read<<cm_all[i]<<"  "<<Yield_all[i]<<"  "<<endl;
	}
    read.close();


    TString ofn = "../simulation.root";
    TFile *fout = new TFile(ofn, "recreate");
    tree->Write();
    h_strip->Write();
    fout->Close();

    cout << endl;
    cout << endl;
    cout << "<Created> ./simulation.root" << endl;
    cout << "...simulation completed!" << endl;
    return 0;
}
