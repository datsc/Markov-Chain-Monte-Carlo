/*****************************************************************************

      program to read in ascii file and write out into binary format

*****************************************************************************/
/*make veryclean make -f make_tipsy*/


#include"tipsydefs.std.h" //struct definitions
#include "xdrfuncs.h" //xdr functions
#include <iostream>
using namespace std;
#include <rpc/rpc.h> 
#include <cmath>

/***************************************************************************

    Helper function

****************************************************************************/


extern "C" void exit(int);

inline void error(const char *s){
    cerr << "error: " << s << '\n';
    exit(1);
}


double write_halo(XDR xdrs,FILE *infile,double eps){

  struct dark_particle dp;
  double m, x,y,z,vx,vy,vz;  

  while (EOF != fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf", &m, &x, &y, &z, &vx, &vy, &vz)) {

    dp.mass = m;
    dp.pos[0] = x;
    dp.pos[1] = y;
    dp.pos[2] = z;
    dp.vel[0] = vx;
    dp.vel[1] = vy;
    dp.vel[2] = vz;
    dp.eps = eps;
    dp.phi = 0.;
    
    xdr_dark(&xdrs, &dp);
  }
  return m;
}


void write_stars(XDR xdrs, FILE *infile,double eps) {

  struct star_particle sp; 
  double m, x,y,z,vx,vy,vz;

   while (EOF != fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf", &m, &x, &y, &z, &vx, &vy, &vz)) {

    sp.mass = m;
    sp.pos[0] = x;
    sp.pos[1] = y;
    sp.pos[2] = z;
    sp.vel[0] = vx;
    sp.vel[1] = vy;
    sp.vel[2] = vz;
    sp.eps =eps;
    sp.phi = 0.;
    
       xdr_star(&xdrs, &sp);
  }
}





/********************************************************************************

     main function

*******************************************************************************/

int main (int argc, char **argv){

  //ifstream ;
  XDR xdrs;
  char *in1, *in2, *out;
  struct dump head;
  double R_h,rb;


  if(argc!=3){
    printf("Forgot the argument:scalelength");
  }else{
    R_h=atof(argv[1]);
    rb=atof(argv[2]);
  }
  

  in1="halo";
  in2="bulge";
  out = "initial.bin";

  FILE *outfile = fopen(out,"w");
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
  if(!outfile) 
    error("The outputfile could not be created!!");

  FILE *infile1= fopen(in1,"r");
  if(!infile1) 
    error("The first file could not be opened!");
  FILE *infile2 = fopen(in2,"r"); 
  if(!infile2) 
    error("The first file could not be opened!");


  //fill header
  head.ndim=3;
  int n_h,n_b;
  float buf;
  fscanf(infile1, "%i %f", &n_h, &buf);
  fscanf(infile2, "%i %f", &n_b, &buf);
  head.time=0.;
  head.ndark=n_h;
  head.nstar=n_b;
  head.nsph=0;
  head.nbodies=head.ndark+head.nstar+head.nsph;
  printf("%s %s %f %f %i %i\n",in1,in2,R_h,rb,n_h,n_b);
  printf("ndark=%i nstar=%i",head.ndark,head.nstar);

   xdr_header(&xdrs, &head);
  
  double eps_h=0.05;
  double m_h; 
  printf("softening halo=%f\n",eps_h);
  m_h = write_halo(xdrs,infile1,eps_h);

  double eps_b=0.05;
  write_stars(xdrs,infile2,eps_b);
  printf("softening bulge=%f\n",eps_b);

  return 0;

}


  
