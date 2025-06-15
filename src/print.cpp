#include <iostream>
#include <cmath>
#include <complex>
#include <time.h>
#include <fstream>
#include <iomanip>


static double rms,dxs,dys,drs,di, st_norm;
static mat vortex_vel(N,2,fill::zeros);

void printVortex(std::vector<Vortex> vortex, double runTime, int fileNumber){


		ostringstream out_vortex;
    	out_vortex << "./data/vortex_xy." << setw(5) << setfill('0') << fileNumber << ends;                  //creates file name for outputting da    ta at time slice
    	string filename_vortex = out_vortex.str();
    	ofstream fout_vortex(filename_vortex.c_str());
    	fout_vortex.precision(12);
 		fout_vortex << scientific;
    	for(int i =0; i < N; i++){
		   fout_vortex << vortex[i].x << " " << vortex[i].y << " " << vortex[i].circ << " " << vortex[i].u << " " << vortex[i].v << endl;
   		}
     
     	fout_vortex.close();

     	ofstream fout_count("./data/curframe.dat");
        fout_count << scientific;
        fout_count.precision(12);
        fout_count << runTime << " " << fileNumber << endl;

	return;
}


void printDiagnostics(double hamiltonianValue, double momentumxValue, double momentumyValue,double momentumangularValue, double runTime, int fileNumber){

    //record values of the Hamiltonian
    ostringstream out_hamiltonian;
    out_hamiltonian << "./output/hamiltonian." << setw(5) << setfill('0') << fileNumber << ends;  //creates file name for outputting data at time slice
    string filename_hamiltonian = out_hamiltonian.str();
    ofstream fout_hamiltonian(filename_hamiltonian.c_str());
    fout_hamiltonian << scientific;
    fout_hamiltonian.precision(12);

    fout_hamiltonian << runTime << " " << hamiltonianValue << " " << momentumxValue << " " << momentumyValue << " " << momentumangularValue << endl;

    fout_hamiltonian.close();

return;

}



void printParameters(){
ofstream fout_data("./parameters.txt");
fout_data.precision(12);
fout_data << "Np = " << Np << " Number of positive vortices" <<  endl;
fout_data << "Nm = " << Nm << " Number of negative vortices" << endl;
fout_data << "N =  " << N << " Total number of vortices" << endl;
fout_data << "Lx = " << Lx << endl;
fout_data << "Ly = " << Ly << endl;
fout_data << "L = " << sqrt(Lx*Ly/(double) N) << " Mean intervortex spacing " << endl;
fout_data << "=========== time stepping =========" << endl;
fout_data << "dtMax = " << dtMax << endl;
fout_data << "rkTolerance = " << rkTolerance << endl;
fout_data << "totalSteps =  " << totalSteps << endl;
fout_data << "=========== boundary =========" << endl;
fout_data << "FLAG_BOUNDARY = " << FLAG_BOUNDARY << endl;
fout_data << "=========== output =========" << endl;
fout_data << "FLAG_OUTPUT = " << FLAG_OUTPUT << endl;
fout_data << "OutputTime = " << outputTime << endl;

//fout_data << "FLAG_STRUCTURE = " << FLAG_STRUCTURE << endl;

fout_data << "=========== initial condition =========" << endl;
fout_data << "FLAG_INITIAL_CONDITION = " << FLAG_INITIAL_CONDITION << endl;
fout_data << "=========== backgorund flow ===========" << endl;
fout_data << "FLAG_BACKGROUND_FLOW = " << FLAG_BACKGROUND_FLOW << endl;
fout_data << "backgroundFlowStrength =  " << backgroundFlowStrength << endl;
fout_data << "=========== add / remove ===========" << endl;
fout_data << "FLAG_VORTEX_REMOVE = " << FLAG_VORTEX_REMOVE << endl;
fout_data << "removalDistanceLower = " << removalDistanceLower << endl;
fout_data << "removalDistanceUpper = " << removalDistanceUpper << endl;
fout_data << "probabilityRemoval = " << probabilityRemoval << endl;
fout_data << "FLAG_VORTEX_ADD = " << FLAG_VORTEX_ADD << endl;
fout_data << "============parallelism============" << endl;
fout_data << "number of threads = " << num_threads << endl;

//output data in the terminal
cout << "Np = " << Np << " Number of positive vortices" <<  endl;
cout << "Nm = " << Nm << " Number of negative vortices" << endl;
cout << "N =  " << N << " Total number of vortices" << endl;
cout << "Lx = " << Lx << endl;
cout << "Ly = " << Ly << endl;
cout << "L = " << sqrt(Lx*Ly/(double) N) << " Mean intervortex spacing " << endl;
cout << "=========== time stepping =========" << endl;
cout << "dtMax = " << dtMax << endl;
cout << "rkTolerance = " << rkTolerance << endl;
cout << "totalSteps =  " << totalSteps << endl;
cout << "=========== boundary =========" << endl;
cout << "FLAG_BOUNDARY = " << FLAG_BOUNDARY << endl;
cout << "=========== output =========" << endl;
cout << "FLAG_OUTPUT = " << FLAG_OUTPUT << endl;
cout << "OutputTime = " << outputTime << endl;

//fout_data << "FLAG_STRUCTURE = " << FLAG_STRUCTURE << endl;

cout << "=========== initial condition =========" << endl;
cout << "FLAG_INITIAL_CONDITION = " << FLAG_INITIAL_CONDITION << endl;
cout << "=========== backgorund flow ===========" << endl;
cout << "FLAG_BACKGROUND_FLOW = " << FLAG_BACKGROUND_FLOW << endl;
cout << "backgroundFlowStrength =  " << backgroundFlowStrength << endl;
cout << "=========== add / remove ===========" << endl;
cout << "FLAG_VORTEX_REMOVE = " << FLAG_VORTEX_REMOVE << endl;
cout << "removalDistanceLower = " << removalDistanceLower << endl;
cout << "removalDistanceUpper = " << removalDistanceUpper << endl;
cout << "probabilityRemoval = " << probabilityRemoval << endl;
cout << "FLAG_VORTEX_ADD = " << FLAG_VORTEX_ADD << endl;
cout << "============parallelism============" << endl;
cout << "number of threads = " << num_threads << endl;

return;
}




/*
void output_vel(mat vortex_vel, double runtime, int counter){


//===============  compute rms of the velocity ===========================
rms = 0.0;

for(int i=0; i < N; i++){

    rms += pow( pow(vortex_vel(i,0),2.0) + pow(vortex_vel(i,1),2.0) , 0.5);

}
rms /= double(N);


ostringstream out_vel;
out_vel << "./Output/rms_velocity." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
string filename_vel = out_vel.str();
ofstream fout_vel(filename_vel.c_str());
fout_vel << scientific;
fout_vel.precision(12);

fout_vel << runtime << " " << rms << endl;
      
fout_vel.close();


  return;
}


void output_structure(mat vortex_xy, rowvec vortex_G, mat structure_function, int counter, irowvec structure_counter){


vortex_vel = step_vortices(vortex_xy, vortex_G);
st_norm=0.0;
for(int i =0; i < N; i++){
	for(int j =i+1; j < N; j++){

		dxs = vortex_xy(j,0) - vortex_xy(i,0);
		dys = vortex_xy(j,1) - vortex_xy(i,1);

    if(dxs >= Lx/2) dxs -= Lx;
    else if(dxs < -Lx/2) dxs += Lx;
    if(dys >= Lx/2) dys -= Ly;
    else if(dys < -Ly/2) dys += Ly;

		drs = pow( dxs*dxs + dys*dys  ,0.5);


		di = int(drs*double(NBIN)/Lx);
		structure_function(di,0) = pow( (vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys,3.0);
		structure_function(di,1) = ((vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys)*pow(vortex_G(j)-vortex_G(i),2.0);
		structure_function(di,2) += pow( (vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys,3.0);
		structure_function(di,3) += ((vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys)*pow(vortex_G(j)-vortex_G(i),2.0);
		structure_counter(di) += 1;
	}
}


ostringstream out_struct;
out_struct << "./Output/structure." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
string filename_struct = out_struct.str();
ofstream fout_struct(filename_struct.c_str());
fout_struct << scientific;
fout_struct.precision(12);

for(int i =0; i < NBIN; i++){
	if(structure_counter(i) < 1){
		st_norm=1.0;
	}
	else{
	 st_norm = double(structure_counter(i));
	}
	fout_struct << i*Lx/double(NBIN) << " " << structure_function(i,0) << " " << structure_function(i,1) << " " << structure_function(i,2)/st_norm << " " << structure_function(i,3)/st_norm << endl;
}
fout_struct.close();


  return;
}
*/

void printRemoval(double runTime, int fileNumber,int & removeLowerCounter, int & removeUpperCounter, double &energyLowerRemoval, double & energyUpperRemoval,double & energyInjection){

      ostringstream out_remove;
      out_remove << "./output/remove." << setw(5) << setfill('0') << fileNumber << ends;                  //creates file name for outputting data at time slice
      string filename_remove = out_remove.str();
      ofstream fout_remove(filename_remove.c_str());
      fout_remove << scientific;
      fout_remove.precision(12);

      fout_remove << runTime << " " << removeLowerCounter << " " << removeUpperCounter << " " << energyLowerRemoval << " " << energyUpperRemoval << " " << energyInjection << endl;
       fout_remove.close();
      
      removeLowerCounter = 0;
      removeUpperCounter = 0;
      energyLowerRemoval = 0.0;
      energyUpperRemoval = 0.0;
      energyInjection = 0.0;

return;


}

