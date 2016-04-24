/* Position angle code and writing our coordinates along amjor axis.*/
#include <stdio.h>
#include <math.h>
#include "/home/dc-ural1/Num_Recipes/nrutil.h"
#include "/home/dc-ural1/Num_Recipes/nrutil.c"

double pi;
double *e1_hat, *e2_hat, *e3_hat;
double *x_helio, *y_helio, *z_helio;
double R_sun;

void find_pos_angle();
void find_xyz_helio();
void find_e1e2e3();
void set_up_xyz_helio_axes();
void set_up_sky_coords();
void find_XYZ_gal();
void find_maj_axis_position();
void find_vert_e2_axis_position();

main(int argc, char *argv[])
{ 
    FILE *fin,*fout,*fcomf; 

  char outfile[100], infile[100],fcom[100];
  double l_cen, b_cen;
  double *X,*Y,*Z,*buf,*vx,*vy,*vz,*eps,*ms,*vlos_star,*rstar;
  double Xcen,Ycen,Zcen,vXcen,vYcen,vZcen;
  double l,b,dlos,xveclos,yveclos,zveclos;
  double *x_star, *y_star, *z_star, x_cen,y_cen,z_cen;
  double *e1_star, *e2_star, *e3_star;
  double e1_cen, e2_cen, e3_cen;
  double *PA, PA_maj;
  double *maj_pos_star,*min_pos_star;
  int N,i;


  sscanf(argv[1],"%d",&N); 
  sscanf(argv[2],"%s",&infile);
  sscanf(argv[3],"%s",&outfile);
  sscanf(argv[4],"%s",&fcom);
   ms=dvector(1,N);X=dvector(1,N);Y=dvector(1,N);Z=dvector(1,N);vx=dvector(1,N);vy=dvector(1,N);vz=dvector(1,N);eps=dvector(1,N);vlos_star=dvector(1,N);rstar=dvector(1,N);
  PA=dvector(1,N);maj_pos_star=dvector(1,N);min_pos_star=dvector(1,N);x_star=dvector(1,N);y_star=dvector(1,N);z_star=dvector(1,N);e1_star=dvector(1,N);e2_star=dvector(1,N);e3_star=dvector(1,N);

  pi = acos(-1.0);
  R_sun = 8.5;

  l_cen = 260.113182*pi/180.0;
  b_cen = -22.222879*pi/180.0;
  
 fcomf = fopen(fcom,"r");
  
fscanf(fcomf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&Xcen,&Ycen,&Zcen,&vXcen,&vYcen,&vZcen,&dlos,&xveclos,&yveclos,&zveclos); 
  fclose(fcomf);
 
 PA_maj = 65.0*pi/180.0; 

  set_up_xyz_helio_axes();
 
   find_xyz_helio(Xcen, Ycen, Zcen, &x_cen, &y_cen, &z_cen);
   set_up_sky_coords(x_cen,y_cen,z_cen); 
   find_e1e2e3(x_cen, y_cen, z_cen, &e1_cen, &e2_cen, &e3_cen);
 
   fin = fopen(infile,"r");
 
     for(i=1;i<=N;i++)
     { 
	 fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf ",&ms[i],&X[i],&Y[i],&Z[i],&vx[i],&vy[i],&vz[i],&eps[i]);
	 vlos_star[i]=((vx[i]-vXcen)*xveclos+(vy[i]-vYcen)*yveclos+(vz[i]-vZcen)*zveclos); 
     }
  
  fclose(fin);
  i=1;

  for(i=1;i<=N;i++){
 find_xyz_helio(X[i], Y[i], Z[i], &x_star[i], &y_star[i], &z_star[i]);
  find_e1e2e3(x_star[i], y_star[i], z_star[i], &e1_star[i], &e2_star[i], &e3_star[i]);
  e1_star[i] -= e1_cen;
  e2_star[i] -= e2_cen;
  e3_star[i] -= e3_cen;
  find_maj_axis_position(e1_star[i], e2_star[i], e3_star[i],PA_maj,&maj_pos_star[i]);
   find_vert_e2_axis_position(e1_star[i], e2_star[i], e3_star[i],PA_maj,&min_pos_star[i]);
   find_pos_angle(e1_star[i], e2_star[i], e3_star[i], &PA[i]);
}
  fout = fopen(outfile,"w");
     for(i=1;i<=N;i++)
     { 
  	 fprintf(fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",ms[i],X[i],Y[i],Z[i],vx[i],vy[i],vz[i],eps[i],maj_pos_star[i],maj_pos_star[i]/101.0*206265.0/60.0,min_pos_star[i],vlos_star[i]);
     }
  fclose(fout);
}
void find_maj_axis_position(double e1, double e2, double e3, double PA, double *pos)
{
  *pos = e2*sin(PA) + e3*cos(PA);
  return;
}
void find_vert_e2_axis_position(double e1, double e2, double e3, double PA, double *pos2)
{
    *pos2 = e2*cos(PA) - e3*sin(PA);

  return;
}
void find_pos_angle(double e1, double e2, double e3, double *PA)
{
  double sin_PA, cos_PA;

  cos_PA = e1*e3_hat[1] + e2*e3_hat[2] + e3*e3_hat[3];
  sin_PA = e1*e2_hat[1] + e2*e2_hat[2] + e3*e2_hat[3];
  *PA = atan2(sin_PA,cos_PA);

  return;
}

void find_XYZ_gal(double l, double b, double d_los, double *X, double *Y, double *Z)
{
  *X = R_sun - d_los*cos(b)*cos(l);
  *Y = -1.0*d_los*cos(b)*sin(l);
  *Z = d_los*sin(b);
 
  return;
}

void find_xyz_helio(double X, double Y, double Z, double *x_star, double *y_star, double *z_star)
{
  *x_star = X*x_helio[1] + Y*x_helio[2] + Z*x_helio[3];
  *y_star = X*y_helio[1] + Y*y_helio[2] + Z*y_helio[3];
  *z_star = X*z_helio[1] + Y*z_helio[2] + Z*z_helio[3];

  return;
 }

void find_e1e2e3(double x_star, double y_star, double z_star, double *e1_star, double *e2_star, double *e3_star)
 {
  *e1_star = x_star*e1_hat[1] + y_star*e1_hat[2] + z_star*e1_hat[3];
  *e2_star = x_star*e2_hat[1] + y_star*e2_hat[2] + z_star*e2_hat[3];
  *e3_star = x_star*e3_hat[1] + y_star*e3_hat[2] + z_star*e3_hat[3];

   return;
}

 void set_up_xyz_helio_axes()
{
  double l_x, l_y, l_z;
  double b_x, b_y, b_z;
  double x_norm, y_norm, z_norm;
 
  /* Set l,b values for heliocentric x,y,z from book */
  /* x axis has RA=0, Dec=0 */
  l_x = 96.337270*pi/180.0;
  b_x = -60.188552*pi/180.0;

  /* y axis has RA=90, Dec=0 */
  l_y = 206.989123*pi/180.0;
  b_y = -11.424494*pi/180.0;

  /* z axis has  Dec=90 */
  l_z = 122.931918*pi/180.0;
  b_z = 27.128251*pi/180.0;

  /* Now define unit vectors in these directions */
  x_helio = dvector(1,3);
  y_helio = dvector(1,3);
  z_helio = dvector(1,3);

  x_helio[1] = R_sun - 1.0E10*cos(b_x)*cos(l_x);
  x_helio[2] = -1.0*1.0E10*cos(b_x)*sin(l_x);
  x_helio[3] = 1.0E10*sin(b_x);
  x_norm = sqrt(x_helio[1]*x_helio[1] + x_helio[2]*x_helio[2] + x_helio[3]*x_helio[3]);
  x_helio[1] /= x_norm;
  x_helio[2] /= x_norm;
  x_helio[3] /= x_norm;

  y_helio[1] = R_sun - 1.0E10*cos(b_y)*cos(l_y);
  y_helio[2] = -1.0*1.0E10*cos(b_y)*sin(l_y);
  y_helio[3] = 1.0E10*sin(b_y);
  y_norm = sqrt(y_helio[1]*y_helio[1] + y_helio[2]*y_helio[2] + y_helio[3]*y_helio[3]);
  y_helio[1] /= y_norm;
  y_helio[2] /= y_norm;
  y_helio[3] /= y_norm;

  z_helio[1] = R_sun - 1.0E10*cos(b_z)*cos(l_z);
  z_helio[2] = -1.0*1.0E10*cos(b_z)*sin(l_z);
  z_helio[3] = 1.0E10*sin(b_z);
  z_norm = sqrt(z_helio[1]*z_helio[1] + z_helio[2]*z_helio[2] + z_helio[3]*z_helio[3]);
  z_helio[1] /= z_norm;
  z_helio[2] /= z_norm;
  z_helio[3] /= z_norm;

  return;
}

void set_up_sky_coords(double x, double y, double z)
{
  double r, R;

  r = sqrt(x*x + y*y + z*z);
  R = sqrt(x*x + y*y);

  e1_hat = dvector(1,3);
  e1_hat[1] = x/r;
  e1_hat[2] = y/r;
  e1_hat[3] = z/r;
  
  e2_hat = dvector(1,3);
  e2_hat[1] = -y/R;
  e2_hat[2] = x/R;
  e2_hat[3] = 0.0;

  e3_hat = dvector(1,3);
  e3_hat[1] = -x*z/(r*R);
  e3_hat[2] = -y*z/(r*R);
  e3_hat[3] = (x*x+y*y)/(r*R);

  return;
}



