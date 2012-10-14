#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multiroots.h>
#include "definitions.h"

int main (int argc, char **argv) {
	
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;
  const size_t n = 7;
  double nbar_max=1., dnbar=0.001 ; 
  double press=0, rho=0 ; 
  
  nbar = atof(argv[1]) ; 

  gsl_multiroot_function f = {&system, n};

  // intial values   
  double x_init[n] = {0.006, 
					  0.005,
					 -0.001,
					  0.85,
					  0.0001,
					 -0.001, 
					 pow(nbar*3.*M_PI*M_PI/(1.
					 + pow(x_init[4], 3.)/n_unit), 
					 1./3.)}; 
					 	  
  gsl_vector *x = gsl_vector_alloc (n);
  
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);  
  gsl_vector_set (x, 2, x_init[2]);
  gsl_vector_set (x, 3, x_init[3]);
  gsl_vector_set (x, 4, x_init[4]);
  gsl_vector_set (x, 5, x_init[5]);    
  gsl_vector_set (x, 6, x_init[6]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);

  // EOS is printed to stdin   
  printf("#nbar			rho			press\n");

  ofstream fout("species.dat", ios::trunc);
  fout << "#nbar		np		nn		ne		nmu		nFL		nFSO		nFSM		nFSP		nFXO		nFXM" ; 
  fout << "		Mp		Mn	MoL		MoSO		MoSM		MoSP		MoXO		MoXM" << endl;
  fout.close() ; 	
	  	
  while(nbar < nbar_max) { 
		   
	do {
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);
		if(status) break ;		// check if solver is stuck 
		status = gsl_multiroot_test_residual (s->f, 1e-10);
		
	} while (status == GSL_CONTINUE && iter < 1000);

    // use values from the previous step in the next one     
    gsl_vector_set (x, 0, gsl_vector_get (s->x, 0)) ; 
    gsl_vector_set (x, 1, gsl_vector_get (s->x, 1)) ; 
    gsl_vector_set (x, 2, gsl_vector_get (s->x, 2)) ; 
    gsl_vector_set (x, 3, gsl_vector_get (s->x, 3)) ; 
    
    if(fabs(gsl_vector_get (s->x, 4))<SMALLVAL) 
		gsl_vector_set (x, 4, x_init[4]) ;
	else       
		gsl_vector_set (x, 4, gsl_vector_get (s->x, 4)) ; 

    if(fabs(gsl_vector_get (s->x, 5))<SMALLVAL) 
		gsl_vector_set (x, 5, x_init[5]) ;
	else
		gsl_vector_set (x, 5, gsl_vector_get (s->x, 5)) ; 
		
    gsl_vector_set (x, 6, gsl_vector_get (s->x, 6)) ; 

	// equation of state 
	eos(x, rho, press) ; 
	printf("%.12le %.12le %.12le\n", nbar, rho, press) ;  

    iter = 0 ; 	
    nbar += dnbar ;                  

	gsl_multiroot_fsolver_set (s, &f, x);
     
  } 

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return 0;

}

//----------------------------------------------------------------------
// Equation of state 
//----------------------------------------------------------------------

void eos(const gsl_vector *x, double &rho, double &press) {

	double s = gsl_vector_get (x, 0);
	double v = gsl_vector_get (x, 1);
	double r = gsl_vector_get (x, 2);
	double a = gsl_vector_get (x, 3);
	double b = gsl_vector_get (x, 4);
	double f = gsl_vector_get (x, 5);    
	double z = gsl_vector_get (x, 6);
			
    // effective masses      		         
    double fdn  = dn-gs*s;
    double fdp  = dp-gs*s;
    double fdL  = MoL(s,b);
    double fdSO = MoSO(s,b);
    double fdSM = MoSM(s,b);
    double fdSP = MoSP(s,b);
    double fdXO = MoXO(s,b);
    double fdXM = MoXM(s,b);

    double mun = sqrt(z*z + fdn*fdn) + gv*v - 0.5*gr*r;
    double mup = sqrt(pow(a*z, 2) + fdp*fdp) + gv*v + 0.5*gr*r;

    // leptons        
    double mme = mue(s,v,r,a,z) ;
    double kFe  = kF(mme, de, 0.);
    double kFmu = kF(mme, dmu, 0.);
            
    // hyperons         
    double muL = mun;            
    double kFL = kF(muL, MoL(s,b), xv*gv*v+xf*gv*f);
    double nFL = nF(muL, MoL(s,b), xv*gv*v+xf*gv*f);

    double muS  = muL;
    double kFSO = kF(muL, MoSO(s,b), xv*gv*v+xf*gv*f);
    double nFSO = nF(muL, MoSO(s,b), xv*gv*v+xf*gv*f);

    double muSM = 2.0*mun - mup;
    double kFSM =kF(2.0*mun-mup, MoSM(s,b), xv*gv*v-gr*r+xf*gv*f);
    double nFSM =nF(2.0*mun-mup, MoSM(s,b), xv*gv*v-gr*r+xf*gv*f);

    double muSP = mup;
    double kFSP = kF(mup, MoSP(s,b), xv*gv*v+gr*r+xf*gv*f);
    double nFSP = nF(mup, MoSP(s,b), xv*gv*v+gr*r+xf*gv*f);
    
    double kFXO = kF(mun, MoXO(s,b), xvx*gv*v+0.5*gr*r+xfx*gv*f);
    double nFXO = nF(mun, MoXO(s,b), xvx*gv*v+0.5*gr*r+xfx*gv*f);

    double kFXM = kF(2.0*mun-mup, MoXM(s,b), xvx*gv*v-0.5*gr*r+xfx*gv*f);
    double nFXM = nF(2.0*mun-mup, MoXM(s,b), xvx*gv*v-0.5*gr*r+xfx*gv*f);
    
    // Number densities: 
    // protons & neutrons 
    double nn  = n(z)*n_unit ;
    double np  = pow(a, 3.)*nn ; 
    // leptons 
    double ne  = nF(mme, de, 0.)*n_unit ; 
    double nmu = nF(mme, dmu, 0.)*n_unit ;  
        
    // equation of state         
    rho = chi0(fdn, z) + chi0(fdp, a*z) + chi0(de, kFe) + chi0(dmu, kFmu) 
        + chi0(fdL, kFL) + chi0(fdXM, kFXM) + chi0(fdXO, kFXO) 
        + chi0(fdSM, kFSM) + chi0(fdSO, kFSO) + chi0(fdSP, kFSP)
        + 0.5*pow(momega/M, 2)*v*v
        + 0.5*pow(mrho/M, 2)*r*r + uo(s)
        + 0.5*pow(mb/M, 2)*b*b + 0.5*pow(mf/M, 2)*f*f
        + 3.*Lv*pow(gr*gv*r*v, 2.) + 0.75*c3*(pow(v, 4) + pow(r, 4));
        
    press = phi0(fdn,z) + phi0(fdp, a*z) + phi0(de, kFe) + phi0(dmu, kFmu)
          + phi0(fdXM, kFXM) + phi0(fdXO, kFXO) + phi0(fdL, kFL) 
          + phi0(fdSM, kFSM) + phi0(fdSO, kFSO) + phi0(fdSP, kFSP)
          + 0.5*pow(momega/M, 2)*v*v
          + 0.5*pow(mrho/M, 2)*r*r - uo(s)
          - 0.5*pow(mb/M, 2)*b*b + 0.5*pow(mf/M, 2)*f*f
          + Lv*pow(gr*gv*r*v, 2.) + 0.25*c3*(pow(v, 4) + pow(r, 4));
	
	 // species output
	 ofstream fout("species.dat", ios::app);	  
     fout.precision(8); 
     fout.scientific; 
     fout.width(14); 
     fout << nbar 			<< " "; fout.width(14);  
     fout << np   			<< " "; fout.width(14); 
     fout << nn   			<< " "; fout.width(14);
     fout << ne   			<< " "; fout.width(14); 
     fout << nmu			<< " "; fout.width(14);
     fout << nFL*n_unit		<< " "; fout.width(14);
     fout << nFSO*n_unit	<< " "; fout.width(14);
     fout << nFSM*n_unit	<< " "; fout.width(14);
     fout << nFSP*n_unit	<< " "; fout.width(14);
     fout << nFXO*n_unit	<< " "; fout.width(14);
     fout << nFXM*n_unit	<< " "; fout.width(14);
     fout << fdp			<< " "; fout.width(14);
     fout << fdn			<< " "; fout.width(14);       
     fout << fdL			<< " "; fout.width(14); 
     fout << fdSO			<< " "; fout.width(14); 
     fout << fdSM			<< " "; fout.width(14); 
     fout << fdSP			<< " "; fout.width(14); 
     fout << fdXO			<< " "; fout.width(14); 
     fout << fdXM 			<< endl;
     fout.close(); 

}

//----------------------------------------------------------------------
// Build the system of equations to solve
//----------------------------------------------------------------------

int system (const gsl_vector *x, 
			void *params, 
			gsl_vector *ff) {

	double s = gsl_vector_get (x, 0);
	double v = gsl_vector_get (x, 1);
	double r = gsl_vector_get (x, 2);
	double a = gsl_vector_get (x, 3);
	double b = gsl_vector_get (x, 4);
	double f = gsl_vector_get (x, 5);    
	double z = gsl_vector_get (x, 6);  


  	double  fdn, fdp, fdL, 
		mun, mup, mme, muL, muS, muSM, muSP,  
		kFL, kFSP, kFSO, nFSO, kFSM, nFSM,  
		nFL, nFSP, kFXO, nFXO, kFXM, nFXM; 
								
    fdn = dn - gs*s;
    fdp = dp - gs*s;
    fdL = MoL(s,b);

    // chemical potentials of neutrons and protons
    mun = sqrt(z*z+fdn*fdn)+gv*v-0.5*gr*r;
    mup = sqrt(pow(a*z, 2)+fdp*fdp)+gv*v+0.5*gr*r;
 
    // chemical potentials of Lambda and Sigmas
    muL = mun;
    muS = muL;
    muSM = 2.0*mun - mup;
    muSP = mup;
            
    // leptons        
    mme = mue(s,v,r,a,z);
            
    // Fermi momenta and corresponding baryon densities         
    kFL=kF(muL,MoL(s,b),xv*gv*v+xf*gv*f);
    nFL=nF(muL,MoL(s,b),xv*gv*v+xf*gv*f);

    kFSO=kF(muS,MoSO(s,b),xv*gv*v+xf*gv*f);
    nFSO=nF(muS,MoSO(s,b),xv*gv*v+xf*gv*f);

    kFSM=kF(muSM,MoSM(s,b),xv*gv*v-gr*r+xf*gv*f);
    nFSM=nF(muSM,MoSM(s,b),xv*gv*v-gr*r+xf*gv*f);

    kFSP=kF(muSP,MoSP(s,b),xv*gv*v+gr*r+xf*gv*f);
    nFSP=nF(muSP,MoSP(s,b),xv*gv*v+gr*r+xf*gv*f);

    kFXO=kF(mun,MoXO(s,b),xvx*gv*v+0.5*gr*r+xfx*gv*f);
    nFXO=nF(mun,MoXO(s,b),xvx*gv*v+0.5*gr*r+xfx*gv*f);

    kFXM=kF(2.0*mun-mup,MoXM(s,b),xvx*gv*v-0.5*gr*r+xfx*gv*f);
    nFXM=nF(2.0*mun-mup,MoXM(s,b),xvx*gv*v-0.5*gr*r+xfx*gv*f);

    //------------------------------------------------------------------
    // left-hand-sides of the equations 
    //------------------------------------------------------------------ 
       
	double y0 = (pow(ms/M, 2))*s+(g3/M)*s*s+g4*s*s*s
	- gs*(1. - gs*s)*(S(1. - gs*s, z)+S(1. - gs*s, a*z))
	- gsL*MoL(s,b)*S(MoL(s,b),kFL)    
	- gsX*MoXM(s,b)*S(MoXM(s,b),kFXM)-gsX*MoXO(s,b)*S(MoXO(s,b),kFXO)
	- gsS*MoSP(s,b)*S(MoSP(s,b),kFSP)-gsS*MoSM(s,b)*S(MoSM(s,b),kFSM) 
	- gsS*MoSO(s,b)*S(MoSO(s,b),kFSO);

	double y1 = v*pow(momega/M, 2)*(1.0+eta1*gs*s+0.5*eta2*gs*gs*s*s)+c3*v*v*v    
	- gv*nB(a,z)
	- gvL*nFL
	- gvX*(nFXM+nFXO)
	- gvS*(nFSP+nFSM+nFSO)
	+ 2.*Lv*pow(gr*gv, 2)*v*r*r;
          
	double y2 = r*pow(mrho/M, 2)*(1.0+gs*etar*s)+c3*r*r*r
	- 0.5*gr*(n(a*z)-n(z))
	- 0.5*gr*(nFXM-nFXO)
	- gr*(nFSP-nFSM)
	+ 2.*Lv*pow(gr*gv, 2)*v*v*r;      

    // electric charge neutrality 
	double y3 = pow((3*M_PI*M_PI)*(nF(mme, de, 0.) 
	+ nF(mme, dmu, 0.) + nFXM + nFSM - nFSP), 1./3.) - z*a;
    	
	double y4 = pow(mb/M, 2)*b-xbL*gs*(fdL*S(fdL,kFL))
	- xbX*gs*(MoXO(s,b)*S(MoXO(s,b),kFXO)
	+ MoXM(s,b)*S(MoXM(s,b),kFXM))
	- gbS*MoSP(s,b)*S(MoSP(s,b),kFSP)
	- gbS*MoSM(s,b)*S(MoSM(s,b),kFSM)
	- gbS*MoSO(s,b)*S(MoSO(s,b),kFSO);

	double y5 = pow(mf/M, 2)*f-xf*gv*nFL
	- xf*gv*(nFSO+nFSP+nFSM);
    
    // baryon number conservation      
	double y6 = n_unit*(nB(a,z) + nFL + nFSO + nFSP + nFSM + nFXO + nFXM) 
	- nbar ; 

	gsl_vector_set (ff, 0, y0);
	gsl_vector_set (ff, 1, y1);
	gsl_vector_set (ff, 2, y2);
	gsl_vector_set (ff, 3, y3);
	gsl_vector_set (ff, 4, y4);
	gsl_vector_set (ff, 5, y5);
	gsl_vector_set (ff, 6, y6);	
	
	return GSL_SUCCESS;
	
}

//----------------------------------------------------------------------
// Auxilliary functions 
//----------------------------------------------------------------------

double uo(double s) { 
	
     return 0.5*(pow(ms/M, 2))*s*s+(g3/(3.0*M))*pow(s, 3) + 0.25*g4*pow(s, 4);

}

double n(double x) { 
    
    return (1.0/(3.0*(M_PI*M_PI)))*pow(x, 3);
}
    
double chi0(double d, double x) { 
    
    double sqrtd2x2 = sqrt(d*d+x*x);

    return (1.0/(8.0*M_PI*M_PI))*(pow(d, 4)*log(d)+x*sqrtd2x2*(d*d+2.0*x*x)
		- pow(d, 4)*log(x + sqrtd2x2));
    
} 

double phi0(double d, double x) { 

    double sqrtd2x2 = sqrt(d*d+x*x);
    
    return -(1./(8.*M_PI*M_PI))*pow(d, 4)*log(d)
		+ (1./(24.*M_PI*M_PI))*(x*sqrtd2x2*(-3.*d*d+2.*x*x)
		+ 3.*pow(d, 4)*log(x + sqrtd2x2));

}

double S(double d, double x) { 

	double sqrtd2x2 = sqrt(d*d+x*x);
	
    return (1.0/(M_PI*M_PI))*(0.25*d*d*log(d*d)
			+0.5*(x*sqrtd2x2 - (d*d)*log(x + sqrtd2x2)));

}    

double nB(double a, double x) { 
	
    return n(x)+n(a*x);
}
    
double mue(double s, double v, double r, double a, double z) { 

    return sqrt(z*z + pow(dn-gs*s, 2)) 
		- sqrt(pow(a*z, 2) + pow(dp-gs*s, 2))-gr*r;

}

double nF(double mu, double d, double r) { 
   
    double qq = -d*d + pow(mu-r, 2);
    
    if(qq < 0.)	return 0.; 
    else 		return (1./(3.*M_PI*M_PI))*pow(qq, 1.5);

} 

double kF(double mu, double d, double r) { 
	
    double qq = -d*d + pow(mu-r, 2);
    
    if (qq < 0.)	return 0.;
    else 			return sqrt(qq);

} 

double MoL(double s, double b) { 
	
    return dLambda-gs*(xsL*s+xbL*b);

}

double MoSO(double s, double b) { 
	
    return dSigmaO-gs*(xsS*s+xbS*b);
   
}

double MoSM(double s, double b) { 
	
    return dSigmaM-gs*(xsS*s+xbS*b);

} 
    
double MoSP(double s, double b) { 
	
	return dSigmaP-gs*(xsS*s+xbS*b);

}

double MoXO(double s, double b) { 

	return dXiO-gs*(xsX*s+xbX*b);
	
}

double MoXM(double s, double b) { 
    
    return dXiM-gs*(xsX*s+xbX*b);
    
}


int print_state (size_t iter, gsl_multiroot_fsolver *s) {
	
  printf ("iter = %3lu %lf %.6le %.6le %.6le %.6le %.6le %.6le %.6le\n",
          iter,
          nbar,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->x, 4), 
          gsl_vector_get (s->x, 5),
          gsl_vector_get (s->x, 6));

}
