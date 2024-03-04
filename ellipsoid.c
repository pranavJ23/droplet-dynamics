// #include <stdio.h>
// #include <math.h>
// int main() {
//   printf("hello world");
//   float D_0=2.5;   //mm
//   float V_0=2;    //m/s
//   float density=998;  //kg/m3
//   float mu_d=1;      //milliPa
//   float ST=73;       //milliN/m
//   float weber=density*V_0*V_0*D_0/ST;
//   float Re=density*V_0*D_0/mu_d;
//   float Oh=mu_d/sqrt(density*ST*D_0);
//   return 0;
// }

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"
#include "contact.h"
#include "two-phase.h"

// scalar f[], * interfaces = {f};


double theta0 = 150.0;

vector h[];
h.t[back] = contact_angle (theta0*pi/180.);
h.r[back] = contact_angle (theta0*pi/180.);

#ifdef BC
u.t[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);
u.t[front] = dirichlet(0);
u.t[back] = dirichlet(0);

u.n[top] = dirichlet(0);
u.n[bottom] = dirichlet(0);
u.n[left]   = dirichlet(0);
u.n[right]  = dirichlet(0);
u.n[front] = dirichlet(0);
u.n[back] = dirichlet(0);

u.r[top] = dirichlet(0);
u.r[bottom] = dirichlet(0);
u.r[left]   = dirichlet(0);
u.r[right]  = dirichlet(0);
u.r[front] = dirichlet(0);
u.r[back] = dirichlet(0);
#endif

#define MAXLEVEL 7

int main()
{
  rho1=750.0;
  rho2=1.0;
  mu2 = 1.0/30.0;
  mu1 = 100.0/30.0;
  //CFL = 0.1;
  // const face vector muc[] = {.1,.1,.1};
  // mu = muc;
  // #if BC
  f.height = h;
  // #endif
  f.sigma = 50.0;
  N = 1 << MAXLEVEL; 
 // origin(-0.5,-0.5,-0.5);
  init_grid(64);
    run();
  
}


// u.t[back] = dirichlet(0);

event init (t = 0)
{
  fraction (f, - (sq(x)/4.0 + sq(y)/4.0 + sq(z)/1.0 - sq(0.1)));
}

event acceleration(i++){
  
  vector av=a;
  if(t<0.05){
  foreach()
  // theta0=10;
  av.z[]=0;
  }
  else{
    foreach()
    av.z[]=0;
    // theta0=170;
    // foreach()
    // av.z[]-=10*(t/100)+0.1;
  }
}

event logfile (i += 10; t <= 1)
{
  scalar kappa[];
  cstats cs = curvature (f, kappa);
  foreach()
    if (f[] <= 1e-3 || f[] >= 1. - 1e-3)
      kappa[] = nodata;
  stats s = statsf (kappa);

  fprintf (fout, "%g %g %g %g %g %g %g %d %d %d %d\n", t, normf(u.x).max,
	   s.min, s.sum/s.volume, s.max, s.stddev, statsf(f).sum,
	   cs.h, cs.f, cs.a, cs.c);
  fflush (fout);
}

// event output (i+= 100; t <= 5){

//   char filename[100];
//   sprintf(filename,"vof%d.dat",i);
//   FILE* f1;
//   f1 = fopen(filename,"w");
//   output_facets (f,f1);
// }

event end (t = end)
{
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4.*statsf(f).sum;
  fprintf (stderr, "%d %.5g %.3g %.5g %.3g %.5g\n",
	   N, R, s.stddev, V, normf(u.x).max, t);
}

event movie (i += 5; t <= 1)
{
    view (fov = 26.6776, quat = {0.474458,0.144142,0.234923,0.836017},
	  tx = -0.0137556, ty = -0.00718937, bg = {1,1,1},
	  width = 758, height = 552);


    draw_vof ("f");
    draw_vof ("f", edges = true);
    cells (n = {0,0,1}, alpha = 0.0, lc = {1,0,0});
    mirror (n = {1,0,0}) {
      draw_vof ("f");
      draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
    }
    mirror (n = {0,1,0}) {
      draw_vof ("f");
      draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
      mirror (n = {1,0,0}) {
      draw_vof ("f");
      draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
      }
    }
    save ("movie_trial_without_sitting.mp4");
  
}

#if TREE
event adapt (i++) {
#if 1
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1}, (double[]){1e-3}, minlevel = 3, maxlevel = MAXLEVEL);
#else
  adapt_wavelet ({f}, (double[]){1e-4}, minlevel = 3, maxlevel = MAXLEVEL);
#endif
}
#endif
