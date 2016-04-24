/* This code uses J2000 factors to convert proper motions to velocities */
/* => ensure that all input is in J2000 */
/* The output goes to convvel.dat */

#include <stdio.h>
#include <math.h>

static void johnson();
static void convert_proper_motions();

void convert_proper_motions(double ra,double dec,double d_los,double mu_alpha,double mu_delta,double v_los, double *V_x, double *V_y, double *V_z)
 

{
  //printf("%f %f",mu_alpha,mu_delta);
  double U,V,W;
  double pi;

  pi = acos(-1.0);

  /* Assumes ra, dec in degrees */
  ra /= 180.0/pi;
  dec /= 180.0/pi;
 
  /* Assumes proper motions are in mas/yr - convert proper motions to arcsec/yr */
  mu_alpha /= 1000.0;
  mu_delta /= 1000.0;
  
  /* Perform conversion - assumes d_los in kpc */
  johnson(ra,dec,1.0/(d_los*1000.0),mu_alpha,mu_delta,v_los,&U,&V,&W);
  
  /* To conform to definition that U is positive away from Galactic centre */
  U *= -1.0;

  /* Correct for solar motion. */
#ifdef DINESCU
  U += -11.0;
  V += 14.0;
  W += 7.5;
#else /* Use "standard values" for solar peculiar velocity */
   U += -9.0;
  V += 12.0;
  W += 7.0;
#endif

  *V_x = U;
  *V_y = -1.0*(V+220.0);
  *V_z = W;

  return;
}

/* Implementation of transformations given in Johnson & Soderblom (1987, AJ, 93, 864) */
void johnson(double ra, double dec, double parallax, double mu_alpha, double mu_delta, double v_hel, double *U, double *V, double *W)
{
  double k=4.74057;
  double pi;
  pi = acos(-1.0);
#ifdef B1950
  double alpha_ngp=192.25/180.0*pi;  
  double delta_ngp=27.4/180.0*pi;
  double theta_0=123.0/180.0*pi;
#else /* Default to J2000 */
  double alpha_ngp=192.8594812*pi/180.0;
  double delta_ngp=27.1282512*pi/180.0;
  double theta_0=122.9319186*pi/180.0;
#endif 

  *U=v_hel*(cos(delta_ngp)*cos(theta_0)*sin(dec) + cos(dec)*sin(ra)*(-1.0*(cos(theta_0)*sin(alpha_ngp)*sin(delta_ngp)) + cos(alpha_ngp)*sin(theta_0)) + cos(ra)*cos(dec)*(-1.0*(cos(alpha_ngp)*cos(theta_0)*sin(delta_ngp)) - sin(alpha_ngp)*sin(theta_0))) + (k*mu_alpha/parallax)*(cos(ra)*(-1.0*(cos(theta_0)*sin(alpha_ngp)*sin(delta_ngp)) + cos(alpha_ngp)*sin(theta_0)) - sin(ra)*(-1.0*(cos(alpha_ngp)*cos(theta_0)*sin(delta_ngp)) - sin(alpha_ngp)*sin(theta_0))) + (k*mu_delta/parallax)*(cos(delta_ngp)*cos(dec)*cos(theta_0) - sin(ra)*sin(dec)*(-1.0*(cos(theta_0)*sin(alpha_ngp)*sin(delta_ngp)) + cos(alpha_ngp)*sin(theta_0)) - cos(ra)*sin(dec)*(-1.0*(cos(alpha_ngp)*cos(theta_0)*sin(delta_ngp)) - sin(alpha_ngp)*sin(theta_0)));

  *V = (k*mu_alpha/parallax)*(-1.0*(sin(ra)*(cos(theta_0)*sin(alpha_ngp) - cos(alpha_ngp)*sin(delta_ngp)*sin(theta_0))) + cos(ra)*(-1.0*(cos(alpha_ngp)*cos(theta_0)) - sin(alpha_ngp)*sin(delta_ngp)*sin(theta_0))) + (k*mu_delta/parallax)*(cos(delta_ngp)*cos(dec)*sin(theta_0) - cos(ra)*sin(dec)*(cos(theta_0)*sin(alpha_ngp) - cos(alpha_ngp)*sin(delta_ngp)*sin(theta_0)) - sin(ra)*sin(dec)*(-1.0*(cos(alpha_ngp)*cos(theta_0)) - sin(alpha_ngp)*sin(delta_ngp)*sin(theta_0))) + v_hel*(cos(delta_ngp)*sin(dec)*sin(theta_0) + cos(ra)*cos(dec)*(cos(theta_0)*sin(alpha_ngp) - cos(alpha_ngp)*sin(delta_ngp)*sin(theta_0)) + cos(dec)*sin(ra)*(-1.0*(cos(alpha_ngp)*cos(theta_0)) - sin(alpha_ngp)*sin(delta_ngp)*sin(theta_0)));

  *W = (k*mu_alpha/parallax)*(cos(delta_ngp)*cos(ra)*sin(alpha_ngp) - cos(alpha_ngp)*cos(delta_ngp)*sin(ra)) + (k*mu_delta/parallax)*(cos(dec)*sin(delta_ngp) - cos(alpha_ngp)*cos(delta_ngp)*cos(ra)*sin(dec) - cos(delta_ngp)*sin(alpha_ngp)*sin(ra)*sin(dec)) + v_hel*(cos(alpha_ngp)*cos(delta_ngp)*cos(ra)*cos(dec) + sin(delta_ngp)*sin(dec) + cos(delta_ngp)*cos(dec)*sin(alpha_ngp)*sin(ra));

  return;
}

int main(int argc, char *argv[])
{
double V_x, V_y, V_z;
 double ra; sscanf(argv[1],"%lf",&ra);
 double dec; sscanf(argv[2],"%lf",&dec);
 double d_los; sscanf(argv[3],"%lf",&d_los);
 double mu_alpha; sscanf(argv[4],"%lf",&mu_alpha);
 double mu_delta; sscanf(argv[5],"%lf",&mu_delta);
 double mu_los; sscanf(argv[6],"%lf",&mu_los);
 double kpc,myrs;
 FILE *convvel;
 double xcon,ycon,zcon,stepcon,timecon;

convvel=fopen("convvel.dat","w");
kpc=3.0857e16;
myrs=3.15576e13;
 
 xcon=24.5433;
 ycon=92.0189; 
 zcon=-38.1567; 
 stepcon=0.0001;
 timecon=6000.0;

 convert_proper_motions(ra,dec,d_los,mu_alpha,mu_delta,mu_los,&V_x,&V_y,&V_z);
 //fprintf(convvel,"%E %E %E\n",V_x,V_y,V_z);
 fprintf(convvel,"%f %f %f %f %f %f %f %f \n",xcon,ycon,zcon,(-1*V_x/kpc*myrs),(-1*V_y/kpc*myrs),(-1*V_z/kpc*myrs),stepcon,timecon);
 fclose(convvel);

}


