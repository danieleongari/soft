/*******************************************************************************
File SOFT.h is a header file for program SOFT.c.
*******************************************************************************/
#define NX 512 //Number of bins
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr


/* Function prototypes ********************************************************/
void init_param();
void init_prop();
void init_wavefn();
void single_step(int step);
void pot_prop();
void kin_prop();
void periodic_bc();
void calc_energy();
void four1(double data[], unsigned long nn, int isign);
void create_psif();
void update_psi();
void store_psip();
void calc_ekin();
void calc_epot();
void calc_kv();
void calc_vk();
void calc_deltaE();
void calc_norm();
void print_energy(int step, FILE *f1);
void print_wavefn(int step, FILE *f2, FILE *f3, FILE *f4);
void print_pot(FILE *f5);

/* Input parameters ***********************************************************/
double LX;       /* Simulation box length */
double DT;       /* Time discretization unit */
double TT;
int NSTEP;       /* Number of simulation steps */
int NECAL;       /* Interval to calculate energies */
int NNCAL;       /* Interval to calculate the norm */
int pot_type;    /* Type of potential choice */
double X0,S0,E0,K0;   /* Center-of-mass, spread, energy (or velocity) of initial wave packet */
double BH,BW,EH,D,b; /* Parameters for the potentials */
double M;        /* Mass */
/* Arrays **********************************************************************
psi[NX+2][2]:    psi[i][0|1] is the real|imaginary part of the wave function
                 on mesh point i
psif[(NX+1)*2]:  psif[2i|2i+1] is the real|imaginary part of the wave function
		 on mesh point i 
wrk[NX+2][2]:    Work array for a wave function
T[NX+2]:	 T[i] is the kinetic energy at mesh point i
t[NX+2][2]:	 t[i][] is the kinetic propagator on i (real|imaginary part)
v[NX+2]:         v[i] is the potential energy at mesh point i
u[NX+2][2]:      u[i][] is the potential propagator on i (real|imaginary part)
*******************************************************************************/
double psi[NX+1][2];
double psip[NX+1][2]; 
double psif[(NX+1)*2+1]; //psif[0] is skipped
double wrk[NX+1][2];
double T[NX+1];
double t[NX+1][2];
double v[NX+1];
double u[NX+1][2];

/* Variables *******************************************************************
dx   = Mesh spacing
ekin = Kinetic energy
epot = Potential energy
etot = Total energy
*******************************************************************************/
double dx;
double norm;
double xave;
double kave;
double ekin,epot,etot;
double v2,k2;
double kv,vk,deltaE,E2;
