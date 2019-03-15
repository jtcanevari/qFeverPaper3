# include <Rcpp.h>
# include <vector>


  using namespace Rcpp;

/// nmath.h --------------------------------------------------------------
/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/* Required by C99 but might be slow */
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif

/* To ensure atanpi, cospi,  sinpi, tanpi are defined */
# ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#  define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
# endif

#include <math.h>
#include <float.h> /* DBL_MIN etc */

#include <Rconfig.h>
#include <Rmath.h>

/* Used internally only */
double  Rf_d1mach(int);
double	Rf_gamma_cody(double);

#include <R_ext/RS.h>

/* possibly needed for debugging */
#include <R_ext/Print.h>

/* moved from dpq.h */
#ifdef HAVE_NEARYINT
# define R_forceint(x)   nearbyint()
#else
# define R_forceint(x)   round(x)
#endif
//R >= 3.1.0: # define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-7)
# define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-7*fmax2(1., fabs(x)))

#ifndef MATHLIB_STANDALONE

#include <R_ext/Error.h>
# define MATHLIB_ERROR(fmt,x)		error(fmt,x);
# define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) warning(fmt,x,x2,x3,x4,x5)

#include <R_ext/Arith.h>
#define ML_POSINF	R_PosInf
#define ML_NEGINF	R_NegInf
#define ML_NAN		R_NaN


void R_CheckUserInterrupt(void);
/* Ei-ji Nakama reported that AIX 5.2 has calloc as a macro and objected
to redefining it.  Tests added for 2.2.1 */
#ifdef calloc
# undef calloc
#endif
#define calloc R_chk_calloc
#ifdef free
# undef free
#endif
#define free R_chk_free

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) gettext (String)
#else
#define _(String) (String)
#endif

#else
/* Mathlib standalone */

#include <stdio.h>
#include <stdlib.h> /* for exit */
#define MATHLIB_ERROR(fmt,x)	{ printf(fmt,x); exit(1); }
#define MATHLIB_WARNING(fmt,x)		printf(fmt,x)
#define MATHLIB_WARNING2(fmt,x,x2)	printf(fmt,x,x2)
#define MATHLIB_WARNING3(fmt,x,x2,x3)	printf(fmt,x,x2,x3)
#define MATHLIB_WARNING4(fmt,x,x2,x3,x4) printf(fmt,x,x2,x3,x4)
#define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) printf(fmt,x,x2,x3,x4,x5)

#define ISNAN(x) (isnan(x)!=0)
// Arith.h defines it
#ifndef R_FINITE
#ifdef HAVE_WORKING_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif
#endif
int R_finite(double);

#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	((-1.0) / 0.0)
#define ML_NAN		(0.0 / 0.0)

#define _(String) String
#endif /* standalone */

#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }

/* For a long time prior to R 2.3.0 ML_ERROR did nothing.
We don't report ME_DOMAIN errors as the callers collect ML_NANs into
a single warning.
*/
#define ML_ERROR(x, s) {                                            \
if(x > ME_DOMAIN) {                                                 \
  char *msg = "";                                                   \
  switch(x) {                                                       \
  case ME_DOMAIN:                                                   \
    msg = _("argument out of domain in '%s'\n");	                   \
    break;                                                          \
  case ME_RANGE:                                                    \
    msg = _("value out of range in '%s'\n");	                       \
    break;                                                          \
  case ME_NOCONV:                                                   \
    msg = _("convergence failed in '%s'\n");	                       \
    break;                                                          \
  case ME_PRECISION:                                                \
    msg = _("full precision may not have been achieved in '%s'\n"); \
    break;                                                          \
  case ME_UNDERFLOW:                                                \
    msg = _("underflow occurred in '%s'\n");	                       \
    break;                                                          \
  }                                                                 \
  MATHLIB_WARNING(msg, s);                                          \
}                                                                   \
}

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

/* Formerly private part of Mathlib.h */

/* always remap internal functions */
#define bd0       	Rf_bd0
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define gammalims	Rf_gammalims
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define stirlerr       	Rf_stirlerr
#define pnchisq_raw   	Rf_pnchisq_raw
#define pgamma_raw   	Rf_pgamma_raw
#define pnbeta_raw   	Rf_pnbeta_raw
#define pnbeta2       	Rf_pnbeta2
#define bratio       	Rf_bratio

/* Chebyshev Series */

int	attribute_hidden chebyshev_init(double*, int, double);
double	attribute_hidden chebyshev_eval(double, const double *, const int);

/* Gamma and Related Functions */

void	attribute_hidden gammalims(double*, double*);
double	attribute_hidden lgammacor(double); /* log(gamma) correction */
double  attribute_hidden stirlerr(double);  /* Stirling expansion "error" */

double	attribute_hidden lfastchoose(double, double);

double  attribute_hidden bd0(double, double);

double  attribute_hidden pnchisq_raw(double, double, double, double, double,
                                     int, Rboolean, Rboolean);
double  attribute_hidden pgamma_raw(double, double, int, int);
double	attribute_hidden pbeta_raw(double, double, double, int, int);
double  attribute_hidden qchisq_appr(double, double, double, int, int, double tol);
LDOUBLE attribute_hidden pnbeta_raw(double, double, double, double, double);
double	attribute_hidden pnbeta2(double, double, double, double, double, int, int);

int	Rf_i1mach(int);

/* From toms708.c */
void attribute_hidden bratio(double a, double b, double x, double y,
                             double *w, double *w1, int *ierr, int log_p);


#endif /* MATHLIB_PRIVATE_H */
/// end nmath.h --------------------------------------------------------------


/// --------------------------------------------------------------------------
//define function for pglogis. Distribution function for the logistic distribution

double plogis_cpp(double x, double location, double scale,int lower_tail, int log_p)
{
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(location) || ISNAN(scale))
    return x + location + scale;
#endif
  if (scale <= 0.0)	ML_ERR_return_NAN;
  
  x = (x - location) / scale;
  if (ISNAN(x))	ML_ERR_return_NAN;
  R_P_bounds_Inf_01(x);
  
  if(log_p) {
    // log(1 / (1 + exp( +- x ))) = -log(1 + exp( +- x))
    return -log1pexp(lower_tail ? -x : x);
  } else {
    return 1 / (1 + exp(lower_tail ? -x : x));
  }
}

double pglogis (double q, double location, double scale, double shape, int lower_tail, int log_p) 
{
  double x = (q - location)/scale;
  double out;
  
  if (log_p) {
    if (lower_tail) {
      out = shape * plogis_cpp(x,0,1, TRUE, FALSE);
    }
    else {
      out = log(1 - pow(plogis_cpp(x,0,1, TRUE, FALSE), shape));
    }
  }
  else {
    if (lower_tail) {
      out = exp(shape * plogis_cpp(x,0,1, TRUE, FALSE));
    }
    else {
      out = 1 - pow(plogis_cpp(x,0,1, TRUE, FALSE), shape);
    }
  }
  return out;
}

/// end pglogis 


/// --------------------------------------------------------------------------
/// define dglogis - density function for logistic dist.
double dglogis(double x, double location, double scale, double shape, int logar)
{
  double out;
  x = (x - location)/scale;
  double pdf = log(shape) - log(scale) - x - (shape + 1) * log(1 + exp(-x));
  if (logar){ 
    out = pdf;
  } else {
    out = exp(pdf);
  }
  return out;
}

/// end dglogis 
/// --------------------------------------------------------------------------
/// --------------------------------------------------------------------------


/// --------------------------------------------------------------------------
/// define modulus function (time into calendar day)
double modulus(int a, int b, int c)
{ 
  int out = a % b;
  while(out < c) {
    out = out + b;
  }
  return(out);
}

/// end modulus

/// --------------------------------------------------------------------------
//parameteres for ff

double tp = 253; //day of parturition
double Tg = 150; //gestation period
double tc = tp - Tg; //day conception
double Tl = tp - 4*7; //from here on if infected goes to IP2
double ta = tp - 25; //average day of abortion

/// --------------------------------------------------------------------------

/// define forcing function Qp(t) 

// create hazard function for logistic distribution fitted to Autumn kiddings 2015
double HF_cpp(double x)
{
  double out = exp(dglogis(x, tp, 4.987, 1.125, TRUE) - pglogis(x,tp, 4.987, 1.125, FALSE,TRUE));
  return out;
}

//select value for HF depending on the value of t
// [[Rcpp::export]]
double Qp(double t) 
{
  double out;    
  if (modulus(t,365,1) > (tp+40) || modulus(t,365,1) < (tp-30)) {
    out = 0;
  } else {
    out = HF_cpp(modulus(t,365,1));
  }
  return out;
}

/// --------------------------------------------------------------------------
/// define forcing function Qc(t)

//create hazard function for a logistic dist fitted to kidding dist & moved 150d back from kidding date
double HF2_cpp(double x) {
  double out = exp(dglogis(x, tc,4.987, 1.125, TRUE) - pglogis(x, tc, 4.987, 1.125,FALSE, TRUE));
  return out;
}

//select value depending on t
// [[Rcpp::export]]
double Qc(double t) {
  double out;    
  if (modulus(t,365,1) > (tc+40) || modulus(t,365,1) < (tc-30)) {
    out = 0;
  } else {
    out = HF2_cpp(modulus(t,365,1));
  }
  return out;
}

/// --------------------------------------------------------------------------
/// define forcing function Qa(t)

double HF3 (double x) {
  double out = exp(dglogis(x, ta, 6.645, 1.068, TRUE) - pglogis(x, ta, 6.645, 1.068, FALSE, TRUE)); //guarda que F and T estan dados vuelta con respecto a la version en R
  return out;
}

// [[Rcpp::export]]
double Qa(double t) {
  double out;
  if (modulus(t,365,1) < (ta-(11.44*3)) || modulus(t,365,1) > (ta+(11.44*3))){
    out = 0;
  } else {
    out = HF3(modulus(t,365,1));
  }
  return out;
}

/// --------------------------------------------------------------------------
/// define forcing function Qin(t) 
// [[Rcpp::export]]
double Qin (double t) {
  double out;
  if (modulus(t,365,1) < tp+60 || modulus(t,365,1) > tp+90) {
    out = 0;
  } else {
    out = 0.214;  
  }
  return out;
}

/// --------------------------------------------------------------------------
/// define forcing function Q4(t) 
// [[Rcpp::export]]
double Q4 (double t) {
  double out;
  if ((modulus(t,365,1) > Tl) and (modulus(t,365,1) < tp + 40)){
    out = 0;
  } else {
    out = 1;
  }
  return out;
}

/// --------------------------------------------------------------------------
/// define Gillespie

// [[Rcpp::export]]
Rcpp::List CaneGillespie(double t_start, double t_end, std::vector<double> n_init, std::vector<double> par)
{
  double t=t_start, u=0.0; 
  
  double phi = par[0], mu = par[1], fI = par[2], fJ = par[3], gamma = par[4], alpha = par[5],
    beta = par[6], p = par[7], epsilon_p = par[8], epsilon_f = par[9], muE = par[10]; 
  
  double rates[33], E=n_init[16];
  
  int SNP=n_init[0], SP=n_init[1], INP1=n_init[2], INP2=n_init[3], INP3=n_init[4], INP4=n_init[5], 
    IP=n_init[6], IP2=n_init[7], JNP=n_init[8], JP=n_init[9], RNP=n_init[10], RP=n_init[11], Y=n_init[12],
    A=n_init[13], K=n_init[14], KI=n_init[15], i=0;
  
  int N=SNP+SP+INP1+INP2+INP3+INP4+IP+IP2+JNP+JP+RNP+RP; //total pop size
  
// vectors to store outputs
//  NumericVector t_vec, E_vec;
  std::vector<double> t_vec, E_vec;
  t_vec.push_back(t);
  E_vec.push_back(E);

//  IntegerVector SNP_vec, SP_vec, INP_vec, IP_vec, RNP_vec, RP_vec, N_vec;
  std::vector<int> SNP_vec, SP_vec, INP1_vec, INP2_vec, INP3_vec, INP4_vec, IP_vec, IP2_vec, JNP_vec,
   JP_vec, RNP_vec, RP_vec, Y_vec, A_vec, K_vec, KI_vec, N_vec;

  SNP_vec.push_back(SNP);
  SP_vec.push_back(SP);
  INP1_vec.push_back(INP1);
  INP2_vec.push_back(INP2);
  INP3_vec.push_back(INP3);
  INP4_vec.push_back(INP4);
  IP_vec.push_back(IP);
  IP2_vec.push_back(IP2);
  JNP_vec.push_back(JNP);
  JP_vec.push_back(JP);
  RNP_vec.push_back(RNP);
  RP_vec.push_back(RP);
  Y_vec.push_back(Y);
  A_vec.push_back(A);
  K_vec.push_back(K);
  KI_vec.push_back(KI);
  N_vec.push_back(N);

  // Initialise Rcpp RNG
  
  RNGScope scope;

  // Start the loop
  
  while(t<t_end){
    
    // 1. Calculate cumulative rates
    
    rates[0] = Y*Qin(t);
    rates[1] = (phi*(SNP+SP)-SP)*Qc(t);
    rates[2] = (phi*(INP1+IP/4)-IP/4)*Qc(t);
    rates[3] = (phi*(INP2+IP/4)-IP/4)*Qc(t);
    rates[4] = (phi*(INP3+IP/4)-IP/4)*Qc(t);
    rates[5] = (phi*(INP4+IP/4)-IP/4)*Qc(t);
    rates[6] = (phi*(RNP+RP)-RP)*Qc(t);
    rates[7] = (phi*(JNP+JP)-JP)*Qc(t);
    rates[8] = SP*Qp(t);
    rates[9] = (1-alpha)*IP*Qp(t);
    rates[10] = alpha*IP*Qp(t);
    rates[11] = (1-alpha)*IP2*Qp(t);
    rates[12] = alpha*IP2*Qp(t);
    rates[13] = JP*Qp(t);
    rates[14] = RP*Qp(t);
    rates[15] = mu*SNP;
    rates[16] = mu*INP1;
    rates[17] = mu*INP2;
    rates[18] = mu*INP3;
    rates[19] = mu*INP4;
    rates[20] = mu*JNP;
    rates[21] = mu*RNP;
    rates[22] = (beta*SNP*E)/N;
    rates[23] = Q4(t)*beta*SP*E/N;
    rates[24] = beta*(1-Q4(t))*SP*E/N;
    rates[25] = gamma*INP1;
    rates[26] = gamma*INP2;
    rates[27] = gamma*INP3;
    rates[28] = (1-p)*gamma*INP4;
    rates[29] =  p*gamma*INP4;
    rates[30] = alpha*fI*IP*Qa(t);
    rates[31] =  (1-alpha)*fI*IP*Qa(t);
    rates[32] = fJ*JP*Qa(t);
    
    for (i = 0; i < 33; i++) {
      if(rates[i]<0){rates[i] = 0;}
      if(i>0){rates[i] = rates[i]+rates[i-1];}   
    }

    // 2. Draw the time to next event
    
    double tstep = R::rexp(1/rates[32]); //should add small number to prevent ceros in rate[32]?
    
    if(tstep < 14){
    
      // 3. Draw the next event
    
      u = R::runif(0,1)*rates[32];
      
      i=0;
    
      while(rates[i]<u) i++;
    
        // 4. Apply event i
    
        switch(i){
    
        case 0: Y--; SNP++; break;
        case 1: SNP--; SP++; break;
        case 2: INP1--; IP++; break;
        case 3: INP2--; IP++; break;
        case 4: INP3--; IP++; break;
        case 5: INP4--; IP++; break;
        case 6: RNP--; RP++; break;
        case 7: JNP--; JP++; break;
        case 8: SP--; SNP++; K++; break;
        case 9: IP--; RNP++; K++; KI++; E += epsilon_p; break;
        case 10: IP--; JNP++; K++; KI++; E += epsilon_p; break;
        case 11: IP2--; RNP++; K++; E += epsilon_p; break;
        case 12: IP2--; JNP++; K++; E += epsilon_p; break;
        case 13: JP--; RNP++; K++; KI++; E += epsilon_p; break;
        case 14: RP--; RNP++; K++; break;
        case 15: SNP--; Y++; break;
        case 16: INP1--; Y++; break;
        case 17: INP2--; Y++; break;
        case 18: INP3--; Y++; break;
        case 19: INP4--; Y++; break;
        case 20: JNP--; Y++; break;
        case 21: RNP--; Y++; break;
        case 22: SNP--; INP1++; break;
        case 23: SP--;  IP++; break;
        case 24: SP--; IP2++; break;
        case 25: INP1--; INP2++; break;
        case 26: INP2--; INP3++; break;
        case 27: INP3--; INP4++; break;
        case 28: INP4--; SNP++; break;
        case 29: INP4--; RNP++; break;
        case 30: IP--; JNP++; A++; E += epsilon_p; break;
        case 31: IP--; RNP++; A++; E += epsilon_p; break;
        case 32: JP--; RNP++; A++; E += epsilon_p; break;        
    
        default: printf("\nError in choice of event.\n");
    
      }
    
    int x[13] = {SNP,SP,INP1,INP2,INP3,INP4,IP,IP2,JNP,JP,RNP,RP,Y}; //check if any compartment is negative
    int cnt = 0;
    int l = sizeof(x)/sizeof(x[0]);
    
    for (i = 0; i < l; i++) {
      if(x[i]<0){cnt++;}   
    }
    
    if (cnt == 0){ //then all good, proceed
    
      //write to vectors
      SNP_vec.push_back(SNP);
      SP_vec.push_back(SP);
      INP1_vec.push_back(INP1);
      INP2_vec.push_back(INP2);
      INP3_vec.push_back(INP3);
      INP4_vec.push_back(INP4);
      IP_vec.push_back(IP);
      IP2_vec.push_back(IP2);
      JNP_vec.push_back(JNP);
      JP_vec.push_back(JP);      
      RNP_vec.push_back(RNP);
      RP_vec.push_back(RP);
      Y_vec.push_back(Y);
      N=SNP+SP+INP1+INP2+INP3+INP4+IP+IP2+JNP+JP+RNP+RP;
      N_vec.push_back(N);
      A_vec.push_back(A);
      K_vec.push_back(K);
      KI_vec.push_back(KI);
      
      t += tstep;
      E -= E*tstep*muE; E += (INP1+INP2+INP3+INP4+IP+IP2+JNP+JP)*epsilon_f*tstep; //environmental decay and contamination by feces
      if(E<0){ //environment can never be negative
        E = 0;
        }
      t_vec.push_back(t);
      E_vec.push_back(E);  
      
      } else { //if any compartment went negative don't accept, wait and recalculate
      
      SNP = SNP_vec.back();
      SP = SP_vec.back();
      INP1 = INP1_vec.back();
      INP2 = INP2_vec.back();
      INP3 = INP3_vec.back();
      INP4 = INP4_vec.back();
      IP = IP_vec.back();
      IP2 = IP2_vec.back();
      JNP = JNP_vec.back();
      JP = JP_vec.back();
      RNP = RNP_vec.back();
      RP = RP_vec.back();
      Y = Y_vec.back();
      N = N_vec.back();
      A = A_vec.back();
      K = K_vec.back();
      KI = KI_vec.back();
      
      SNP_vec.push_back(SNP);
      SP_vec.push_back(SP);
      INP1_vec.push_back(INP1);
      INP2_vec.push_back(INP2);
      INP3_vec.push_back(INP3);
      INP4_vec.push_back(INP4);
      IP_vec.push_back(IP);
      IP2_vec.push_back(IP2);
      JNP_vec.push_back(JNP);
      JP_vec.push_back(JP);      
      RNP_vec.push_back(RNP);
      RP_vec.push_back(RP);
      Y_vec.push_back(Y);
      N=SNP+SP+INP1+INP2+INP3+INP4+IP+IP2+JNP+JP+RNP+RP;
      N_vec.push_back(N);
      A_vec.push_back(A);
      K_vec.push_back(K);
      KI_vec.push_back(KI);
      
      tstep = 0.2;
      t += tstep; //move t
      E -= E*tstep*muE; E += (INP1+INP2+INP3+INP4+IP+IP2+JNP+JP)*epsilon_f*tstep; //environmental decay and contamination by feces
      if(E<0){
        E = 0;
        }
      t_vec.push_back(t); //keep compartments as they were
      E_vec.push_back(E_vec.back());
      }
    
    } else { //if the time step was bigger than limit then again, wait and recalculate
      SNP_vec.push_back(SNP_vec.back());
      SP_vec.push_back(SP_vec.back());
      INP1_vec.push_back(INP1_vec.back());
      INP2_vec.push_back(INP2_vec.back());
      INP3_vec.push_back(INP3_vec.back());
      INP4_vec.push_back(INP4_vec.back());
      IP_vec.push_back(IP_vec.back());
      IP2_vec.push_back(IP2_vec.back());
      JNP_vec.push_back(JNP_vec.back());
      JP_vec.push_back(JP_vec.back());
      RNP_vec.push_back(RNP_vec.back());
      RP_vec.push_back(RP_vec.back());
      Y_vec.push_back(Y_vec.back());
      N_vec.push_back(N_vec.back());
      A_vec.push_back(A_vec.back());
      K_vec.push_back(K_vec.back());
      KI_vec.push_back(KI_vec.back());

      tstep = 0.5;
      t += tstep;
      E -= E*tstep*muE; E += (INP1+INP2+INP3+INP4+IP+IP2+JNP+JP)*epsilon_f*tstep; //environmental decay and contamination by feces
      if(E<0){
        E = 0;
      }
      t_vec.push_back(t);
      E_vec.push_back(E);
    }
  }
  
  List ret;
  ret["time"] = t_vec;
  ret["SNP"] = SNP_vec;
  ret["SP"] = SP_vec;
  ret["INP1"] = INP1_vec;
  ret["INP2"] = INP2_vec;
  ret["INP3"] = INP3_vec;
  ret["INP4"] = INP4_vec;
  ret["IP"] = IP_vec;
  ret["IP2"] = IP2_vec;
  ret["JNP"] = JNP_vec;
  ret["JP"] = JP_vec;
  ret["RNP"] = RNP_vec;
  ret["RP"] = RP_vec;
  ret["Y"] = Y_vec;
  ret["E"] = E_vec;
  ret["A"] = A_vec;
  ret["K"] = K_vec;
  ret["KI"] = KI_vec;
  ret["N"] = N_vec; //temporary
  return ret;
  
}