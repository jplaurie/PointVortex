
const int Np =1;		//number of positive vortices
const int Nm =0;		// number of negative vortices
const int N = Np+Nm; // total number

const int M_max = 5;

const double pi = 3.14159265358979;
const double Lx = 2.0*pi;
const double Ly = 2.0*pi;

const double ham_eps = 1.e-4;//1.e-12;
const double vel_eps = 1.e-4;//1.e-12;

//====time stepping=========================
const double maxdt = 0.1;
const double rktol = 1.e-8;//1.e-6;
const int nsteps = 1.e9;

const bool FLAG_OUTPUT = true;
const int outstep = 1000;

const bool FLAG_PERIODIC_BC = true;
//initial condition===================
const int FLAG_INITIAL_CONDITION = 1;  // 0 = start from file designated in curframe.dat,1= create random distibution, 2 = from specific file stated in initilize.cpp, 3 = 4-vortex random condition 
//======parameters for white noise on the the vortex positions=========================
const bool FLAG_NOISE = false;
const double D = 0.001; // noise amplitude


//========parameters for the removal of vortices=====================
const bool FLAG_REMOVE = false;
const bool FLAG_ADD = false;
//const bool FLAG_REMOVE_OPP_SIGN_ONLY = false;
const double remove_dist_low = 0.25*sqrt(Lx*Ly / double(N) );		// remove opposite signed vortices if closer than this
const double remove_dist_high = 10.0*sqrt(Lx*Ly / double(N) );		// remove opposite signed vortices if further than this
const double L = sqrt(Lx*Ly / double(N) ) ;				// reintroduce vortices randomly but at this length apart
const double prob_remove = 1.0;							//probability of removing large-scale dipoles

//==========for opemmp===========================
const int num_threads = 1;


//===========output FLAGS==================================
const bool FLAG_STRUCTURE = false;
const int NBIN = 256;


const bool FLAG_BACKGROUND_SHEAR = false;
const double shearStrength = 1.0;