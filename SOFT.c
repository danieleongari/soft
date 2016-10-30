/*******************************************************************************
Quantum dynamics (QD) simulation, SOFT algorithm.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SOFT.h"


int main(int argc, char **argv) {
  
  
  int step; /* Simulation loop iteration index */
   
  FILE *f1=fopen("energy.dat","w");
  FILE *f2=fopen("psi.dat","w");
  FILE *f3=fopen("psisq.dat","w");
  FILE *f4=fopen("norm.dat","w");
  FILE *f5=fopen("potential.dat","w");
  printf("  completed      norm         etot     \n") ;
  
  init_param();          // Initialize input parameters 
  init_pot_box();       //  (choice: _box,_mol,_harm,_flat)
  init_prop();          // Initialize the kinetic & potential propagators
  init_wavefn(f2,f3,f4); // Initialize the electron wave function 

  print_pot(f5);
  print_wavefn(0,f2,f3,f4); 

  for (step=1; step<=NSTEP; step++) {
    single_step(step);                  // Main routine: propagator      
    if (step%NECAL==0) {
      calc_energy();
      print_energy(step,f1);
    }
    if (step%NNCAL==0) {
      calc_norm();
      print_wavefn(step,f2,f3,f4);
    }
    if (step%(NSTEP/10)==0){
      int progress=step*100/NSTEP;
     printf("  %3d %%         %6.5f      %8.4f   \n", progress, norm, etot) ; 
    }
  }
  
  return 0;
}
/*------------------------------------------------------------------------------*/
void init_param() {
/*------------------------------------------------------------------------------
  Initializes parameters.
------------------------------------------------------------------------------*/
  /* initialize control parameters */
  LX=50.0;   //Box dimension [au]
 //NX see .h //Number of bins [-]
  TT=100;    //Simulation total time [au]
  DT=0.001;  //Timestep [au]
  
  NECAL=10;    // Every print energy 
  NNCAL=100;  // Every print psi, psisq, norm 

  
  dx    = LX/NX; // Calculate the mesh size [-] 
  NSTEP = TT/DT;  // Number of steps [au] (!Lowest integer)
}

/*----------------------------------------------------------------------------*/
void init_pot_flat() {
  int sx;
  double x;

  X0 = 10.0;    //Initial position of the particle [au]
  M  = 1.0;     //Mass of the particle [au]
  K0 = 3.0;     //Initial velocity [au] 
  S0 = 0.5;     //Width of the gaussian (sigma) [au]

  for (sx=1; sx<=NX; sx++) {
    x = dx*sx;
    v[sx] = 0;   //Potential
  }
}

/*----------------------------------------------------------------------------*/
void init_pot_box() {
  int sx;
  double x, k;

  X0 = 12.0;    //Initial position of the particle [au]
  M  = 1.0;     //Mass of the particle [au]
  K0 = 2.0;     //Initial velocity [au] 
  S0 = 1.0;     //Width of the gaussian (sigma) [au]

  BH=10.0;    /* height of central barrier */
  BW=3.0;    /* width  of central barrier */
  EH=100.0;  /* height of edge barrier    */

  for (sx=1; sx<=NX; sx++) {
    x = dx*sx;
    if (sx==1 || sx==2 || sx==NX-1 || sx==NX) /* Construct the edge potential */
      v[sx] = EH;    
    else if (0.5*(LX-BW)<x && x<0.5*(LX+BW)) /* Construct the barrier potential */
      v[sx] = BH;
    else
      v[sx] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
void init_pot_mol() {
  int sx;
  double x;

  X0 = 15.0;    //Initial position of the particle [au]
  M  = 1.0;     //Mass of the particle [au]
  K0 = 0.0;     //Initial velocity [au] 
  S0 = 1.0;     //Width of the gaussian (sigma) [au]

  b = 0.3;
  D = 1.0;

  for (sx=1; sx<=NX; sx++) {
    x = dx*sx;
    v[sx] = D*(exp(-2*b*(x-10))-2*exp(-b*(x-10)))+D;    //Potential
  }
}

/*----------------------------------------------------------------------------*/
void init_pot_harm() {
  int sx;
  double x;

  X0 = 17.0;    //Initial position of the particle [au]
  M  = 1.0;     //Mass of the particle [au]
  K0 =-0.5;     //Initial velocity [au] 
  S0 = 1.0;     //Width of the gaussian (sigma) [au]

  b = 0.01;    //Harmonic constant [au]
  D = LX/2;     //Equilibrium position [au]

  for (sx=1; sx<=NX; sx++) {
    x = dx*sx;
    v[sx] = b*(x-D)*(x-D);    //Potential
  }
}


/*----------------------------------------------------------------------------*/
void init_prop() {
  int sx;
  double k, x;

  /* Set up kinetic propagators */
  for (sx=0; sx<=NX; sx++) {
    if (sx < NX/2)
      k = 2*M_PI*sx/LX;            
    else
      k = 2*M_PI*(sx-NX)/LX;       
    
    /* kinetic operator */ 
    T[sx] = 0.5*k*k/M;
    /* kinetic propagator */
    t[sx][0] = cos(-DT*T[sx]);
    t[sx][1] = sin(-DT*T[sx]);
  }
  
  /* Set up potential propagator */
  for (sx=1; sx<=NX; sx++) {
    x = dx*sx;
    /* Half-step potential propagator */
    u[sx][0] = cos(-0.5*DT*v[sx]);
    u[sx][1] = sin(-0.5*DT*v[sx]);
  }
}

/*----------------------------------------------------------------------------*/
void init_wavefn() {
/*------------------------------------------------------------------------------
  Initializes the wave function as a traveling Gaussian wave packet.
------------------------------------------------------------------------------*/
  int sx,s;
  double x,gauss,psisq,norm_fac;

  /* Calculate the the wave function value mesh point-by-point */
  for (sx=1; sx<=NX; sx++) {
    x = dx*sx;
    gauss = exp(-(x-X0)*(x-X0)/4.0/(S0*S0));
    psi[sx][0] = gauss*cos(K0*M*(x-X0)); 	//Re 
    psi[sx][1] = gauss*sin(K0*M*(x-X0));        //Im
  }

  /* Normalize the wave function */
  psisq=0.0;
  for (sx=1; sx<=NX; sx++) {
      psisq += psi[sx][0]*psi[sx][0];
      psisq += psi[sx][1]*psi[sx][1];
  }
  psisq *= dx;
  norm_fac = 1.0/sqrt(psisq);
  for (sx=1; sx<=NX; sx++) {
      psi[sx][0] *= norm_fac;
      psi[sx][1] *= norm_fac;
  }
  periodic_bc();

}

/*----------------------------------------------------------------------------*/
void periodic_bc() {
/*------------------------------------------------------------------------------
  Applies the periodic boundary condition to wave function PSI, by copying
  the boundary values to the auxiliary array positions at the other ends.
------------------------------------------------------------------------------*/
  int s;

    psi[0][0] = psi[NX][0];
    psi[0][1] = psi[NX][1];
    psi[NX+1][0] = psi[1][0];
    psi[NX+1][1] = psi[1][1];
}

/*----------------------------------------------------------------------------*/
void single_step(int step) {
/*------------------------------------------------------------------------------
  Propagates the electron wave function for a unit time step, DT.
------------------------------------------------------------------------------*/
  int j;
  
  pot_prop();  /* half step potential propagation */
  create_psif();
  four1(psif, (unsigned long) NX, -1);
  for (j=0; j <= 2*(NX+1); j++) 
    psif[j] /= NX;
  update_psi();
  kin_prop(); /* step kinetic propagation   */
  if (step%NECAL==0)
    calc_ekin();
  store_psip();
  create_psif();
  four1(psif, (unsigned long) NX, 1);
  update_psi();
  pot_prop();  /* half step potential propagation */
  if (step%NECAL==0)		
    calc_epot();	
}

/*----------------------------------------------------------------------------*/
void pot_prop() {
/*------------------------------------------------------------------------------
  Potential propagator for a half time step, DT/2.
------------------------------------------------------------------------------*/
  int sx;
  double wr,wi;
  
  for (sx=1; sx<=NX; sx++) {
    wr=u[sx][0]*psi[sx][0]-u[sx][1]*psi[sx][1];
    wi=u[sx][0]*psi[sx][1]+u[sx][1]*psi[sx][0];
    psi[sx][0]=wr;
    psi[sx][1]=wi;
  }
  periodic_bc();
}

/*----------------------------------------------------------------------------*/
void kin_prop() {
/*------------------------------------------------------------------------------
  Kinetic propagation for t step.
-------------------------------------------------------------------------------*/
  int sx,s;
  double wr,wi;
  
  for (sx=1; sx<=NX; sx++) {
    wr=t[sx][0]*psi[sx][0]-t[sx][1]*psi[sx][1];
    wi=t[sx][0]*psi[sx][1]+t[sx][1]*psi[sx][0];
    psi[sx][0]=wr;
    psi[sx][1]=wi;
  }
  periodic_bc();	
}

/*----------------------------------------------------------------------------*/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
    if (j > i) {
      SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
      SWAP(data[j+1],data[i+1]);
    }
    m=nn;
    while (m >= 2 && j > m) { 
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  mmax=2;
  while (n > mmax) { /* Outer loop executed log2 nn times. */
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
      for (i=m;i<=n;i+=istep) {
	j=i+mmax; /* This is the Danielson-Lanczos formula. */
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

/*----------------------------------------------------------------------------*/
void create_psif() {
/*------------------------------------------------------------------------------
  Create an array for the fourier transform.
-------------------------------------------------------------------------------*/
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    psif[2*sx-1] = psi[sx-1][0];
    psif[2*sx] = psi[sx-1][1];
  }
}
/*----------------------------------------------------------------------------*/
void update_psi() {
/*------------------------------------------------------------------------------
  Update psi with its fourier transform.
-------------------------------------------------------------------------------*/
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    psi[sx-1][0] = psif[2*sx-1];
    psi[sx-1][1] = psif[2*sx];
  }
}

/*----------------------------------------------------------------------------*/
void store_psip() {

  int sx;

  for (sx=1; sx<=NX; sx++) {
	  psip[sx][0]=psi[sx][0];
	  psip[sx][1]=psi[sx][1];
  }
  psip[0][0]=psip[NX][0];
  psip[0][1]=psip[NX][1];
  psip[NX+1][0]=psip[1][0];
  psip[NX+1][1]=psip[1][1];
}

/*----------------------------------------------------------------------------*/
void calc_ekin() {

  int sx;
  double k;

  ekin = 0.0;
  kave = 0.0;
  for (sx=1; sx<=NX; sx++) {
    if (sx < NX/2)
      k = 2*M_PI*sx/LX;              
    else
      k = 2*M_PI*(sx-NX)/LX;         
    ekin += 0.5/M*k*k*(psi[sx][0]*psi[sx][0]+psi[sx][1]*psi[sx][1]);
    kave += k*(psi[sx][0]*psi[sx][0]+psi[sx][1]*psi[sx][1])*dx;
  }
  ekin *= dx;
  ekin *= NX;
}		

/*----------------------------------------------------------------------------*/
void calc_epot() {
  int sx;

  epot = 0.0;
  for (sx=1; sx<=NX; sx++) {
    epot += v[sx]*(psi[sx][0]*psi[sx][0]+psi[sx][1]*psi[sx][1]);
  }
  epot *= dx;
}

/*----------------------------------------------------------------------------*/
void calc_energy() {

  etot = ekin+epot;
}

/*----------------------------------------------------------------------------*/
void print_energy(int step, FILE *f1) {
  
 fprintf(f1, "%8i %15.10f %15.10f %15.10f\n", step,ekin,epot,etot); //energy.dat
  
}

/*----------------------------------------------------------------------------*/
void calc_norm() {
/*------------------------------------------------------------------------------
  Calculate the norm, average position and velocity.                 
-------------------------------------------------------------------------------*/
  
 int sx;
 double psisq,psisq2;
 double psipsq,psipsq2;

 norm=0.0;
 xave=0.0;

 for (sx=0; sx<=NX-1; sx++) {
  psisq  = psi[sx][0]*psi[sx][0]+psi[sx][1]*psi[sx][1];
  psisq2 = psi[sx+1][0]*psi[sx+1][0]+psi[sx+1][1]*psi[sx+1][1];
  norm += dx*((psisq2+psisq)/2.0);
  xave += dx*((psisq2+psisq)/2.0)*(dx*(sx+0.5));
 }  
  norm=sqrt(norm);
}

/*----------------------------------------------------------------------------*/
void print_wavefn(int step, FILE *f2, FILE *f3, FILE *f4) {
/*------------------------------------------------------------------------------
  Print wf, squared wf and norm 
-------------------------------------------------------------------------------*/
  
 int sx;
 double x, psisqsx, psipsqsx;

 if (step == 0) {
 fprintf(f2,"# sx x psi_Re psi_Im psip_Re psip_Im, index=step/NNCAL \n");
 fprintf(f3,"# sx x psi^2 psip^2, index=step/NNCAL \n");
 fprintf(f4,"# step norm xave kave\n"); 
 }
 
 fprintf(f2,"\n"); 
 fprintf(f2,"\n"); 
 fprintf(f3,"\n"); 
 fprintf(f3,"\n");

 for (sx=1; sx<=NX; sx++){
  x = dx*sx;
  fprintf(f2,"%8i %15.10f %15.10f %15.10f %15.10f %15.10f\n",
              sx,x,psi[sx][0],psi[sx][1],psip[sx][0],psip[sx][1]); // print psi.dat
  psisqsx = psi[sx][0]*psi[sx][0]+psi[sx][1]*psi[sx][1];
  psipsqsx = psip[sx][0]*psip[sx][0]+psip[sx][1]*psip[sx][1];
  fprintf(f3,"%8i %15.10f %15.10f %15.10f\n",sx,x,psisqsx,psipsqsx); // psisq.dat
 }
 if (step > 0) {
  fprintf(f4,"%8i %15.10f %15.10f  %15.10f\n",step,norm,xave,kave);  // norm.dat
 }
}

/*-------------------------------------------------------------------------------*/
void print_pot(FILE *f5) {
  
  int sx;
  double x;
 
 fprintf(f5,"#index x pot");

 for (sx=1; sx<=NX; sx++){
  x=dx*sx;
  fprintf(f5,"%8i %15.10f %15.10f\n",sx,x,v[sx]);  //potential.dat
  }
}

