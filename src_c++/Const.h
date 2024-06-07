
//============== Model Parameters ==============
const int Np = 1;      //Number of positive vortices
const int Nm = 1;		//Number of negative vortices
const int N = Np+Nm;    //Compute the total number of vortices

const int imageTotal = 5;  //Number of periodic images used in full periodic model (5 seem to be sufficient)

const double pi = 3.14159265358979;
const double Lx = 2.0*pi;       //System size
const double Ly = 2.0*pi;

const double R = 1.0;            //System radius of disc domain

//============== Time Stepping Parameters ==============
const double dtMax = 0.05;                  //Maximum timestep to be used in adaptive timestepping scheme 
const double rkTolerance = 1.e-6;  //Maximum tolerance allowed between the 4th and 5th order Runge-Kutta solutions
const int totalSteps = 10000000;//1.e9;            //Total number of time steps


//============== Boundary Conditions ==============
// "periodic" = fully periodic, model of Weiss and McWilliams 1991
// "periodic-x" = periodic in x only, open in y. adapted from Aref JFM 1995
// "infinite" = open boundary conditions
// "disc" = free slip disc of size R
const std::string FLAG_BOUNDARY = "periodic";   

//============== Output Diagonostics and State ==============
const bool FLAG_OUTPUT = true;
const double outputTime = 0.1;      //runTime at which diagonistics and state are outputted to file
const bool FLAG_STRUCTURE = false;  //Turn on/off the computation of structure functions of state
const int NBIN = 256;


//============== Set Initial Conditions ==============
// "currentFrame" = start from file designated in curframe.dat
// "random" = create random distibution
//  "vortexState" = 4-vortex random condition 
const std::string FLAG_INITIAL_CONDITION = "random";


//============== Background Flow ==============

// "shear"  = shear flow v=[A*(y - Ly/2),0]
// "sinh"  = sinh flow v=[A(sinh(y-Ly/2)), 0]
// "none" = no additional background flow
const std::string FLAG_BACKGROUND_FLOW = "none";
const double backgroundFlowStrength = 0.001; //Background flow amplitude

////============== Random Noise ==============
// "gaussian" = random white noise
// "none" = no added noise
const std::string FLAG_NOISE = "none";
const double noiseStrength = 0.001; //Noise amplitude

//============== Routines for Vortex Removal ==============
const bool FLAG_VORTEX_REMOVE = false; //Turn off/on vortex removal at each time step
const double L = sqrt( Lx * Ly / double(N) ) ;	//compute the mean intervortex distance for periodic box
const double Lr = sqrt( pi * R * R / double(N)); //compute the mean intervortex distance for disc
const double removalDistanceLower = 0.5*Lr;		//Remove opposite signed vortices if nearest neighbour is closer than this distance
const double removalDistanceUpper = 400.0*Lr;		//Remove opposite signed vortices if nearest neighbour is futher than this distance
const double probabilityRemoval = 1.0;			//Probability used in candidate of removal of vortex dipole


const bool FLAG_VORTEX_ADD = false; //Turn on/off reinjecting dipole vortices if removed. Added in randomly at distance L

const bool FLAG_RECORD_ENERGY = false; //Turn on/off computation of Hamiltonian difference from add/remove vortices


//============== Parameters for OpenMP ==============
const int num_threads = 1;      //Define the number of threads
