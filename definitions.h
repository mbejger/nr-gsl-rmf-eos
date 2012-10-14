using namespace std;

// global constants
namespace { 
const double mn    =   939.0;		// MeV
const double M     =   938.919;
const double mp    =   938.272;
const double mu    =   105.568;
const double dmu   =   mu/M;
const double me    =   0.511;
const double de    =   me/M; 
const double ms    =   508.194;		// MeV
const double m0    =   197.329;
const double mrho  =   763.0;		// MeV
const double momega=   782.501;		// MeV
const double mb    =   975.0;
const double mf    =   1020.0;
const double g3    =   10.431*m0;	// MeV
const double g4    =   -28.885;
const double g5    =   0.0;
const double g6    =   0.0;
const double gs    =   10.217;		// g_{\sigma}
const double gv    =   12.868;		// g_{\omega}
const double gr    =   9.88958;		// g_{\rho}
const double dp    =   mp/M;
const double dn    =   mn/M;
const double c3    =   0.0;
const double wer   =   101445.46643926227;
const double mPhi  =   1020.0;
const double mLambda = 1115.63;
const double dLambda = mLambda/M;
const double mSigmaP = 11189.37;
const double dSigmaP = mSigmaP/M;
const double mSigmaM = 11197.43;
const double dSigmaM = mSigmaM/M;
const double mSigmaO = 11192.55;
const double dSigmaO = mSigmaO/M;
const double mXiM    = 11321.3;
const double dXiM    = mXiM/M;
const double mXiO    = 11314.9;
const double dXiO    = mXiO/M;
const double xsS   =   4.7103/gs;
const double xbS   =   5.59753/gs;
const double xsL   =   6.2693/gs;
const double xbL   =   5.59753/gs;
const double xsX   =   3.24217/gs;
const double xbX   =   11.7658237/gs;
const double xv    =   2.0/3.0;
const double xvx   =   1.0/3.0;
const double xf    =   -sqrt(2.0)/3.0;
const double xfx   =   -2.0*sqrt(2.0)/3.0;
const double xd    =   0.0;
const double xdx   =   0.0;
const double gsL   =   xsL*gs;
const double gsS   =   xsS*gs;
const double gsX   =   xsX*gs;
const double gvL   =   xv*gv;
const double gvS   =   xv*gv;
const double gvX   =   xvx*gv;
const double gbL   =   xbL*gs;
const double gbS   =   xbS*gs;
const double gbX   =   xbX*gs;
const double eta1  =   0.0;
const double eta2  =   0.0;
const double etar  =   0.0;
const double Lv    =   0.015;
const double noj   =   0.00139245;
const double n_unit=   pow(M/m0, 3.); 
}	

// baryon density; the independent parameter 
double nbar ; 

// auxilary functions 
double uo(double); 
double n(double); 
double chi0(double, double); 
double phi0(double, double); 
double S(double, double); 
double nB(double, double); 
double mue(double, double, double, double, double); 
double nF(double, double, double); 
double kF(double, double, double);
// chemical potential of hyperons 
double MoL(double, double); 
double MoSO(double, double); 
double MoSM(double, double); 
double MoSP(double, double); 
double MoXO(double, double); 
double MoXM(double, double); 
// equation of state 
void eos(const gsl_vector*, double&, double&) ;
// system of nonlinear equations
int system (const gsl_vector*, void*, gsl_vector*) ; 
// diagnostics printout
int print_state (size_t, gsl_multiroot_fsolver*) ; 

#define SMALLVAL 1.e-12