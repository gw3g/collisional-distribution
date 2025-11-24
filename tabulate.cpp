#include "inc/macros.hpp"
#include "inc/integrators.hpp"
#include "inc/htl.hpp"

using namespace std;

/*--------------------------------------------------------------------*/
// useful integrands:

double _S_integrand(double X, void *params) {
  double a = ((double *)params)[0];
  double eps = a - (1.-1/X);
  return eval_S_tot(eps,params)/sqr(X);
}

double _ST_integrand(double X, void *params) {
  double a = ((double *)params)[0];
  double eps = a - (1.-1/X);
  return eval_S_T(eps,params)/sqr(X);
}

double _SL_integrand(double X, void *params) {
  double a = ((double *)params)[0];
  double eps = a - (1.-1/X);
  return eval_S_L(eps,params)/sqr(X);
}

double _S_times_eps(double X, void *params) {
  double eps = (X);
  return eps*eval_S_tot(eps,params);
}

double _S_subtr(double X, void *params) {
  double eps = 1.-.1*(1.-1./X);
  //double eps = (X);
  return .1*eps*( eval_S_tot(eps,params) - .5/sqr(eps) )/sqr(X);
}

double _S_sin_integrand(double X, void *params) {
  double eps = X;
  //double eps = - (1.-1/X);
  double eta = ((double *)params)[0];
  //double eps = t/eta;
  //cout << " eta = " << eta << " , eps = " << eps << endl;
  return sin(eta*eps)*eval_S_tot(eps,params);
}

double _ST_cos_integrand(double eps, void *params) {
  //double eps = - (1.-1/X);
  double eta = ((double *)params)[0];
  //double eps = t/eta;
  //cout << " eta = " << eta << " , eps = " << eps << endl;
  //return (1.-cos(eta*eps))*eval_S_T(eps,params);
  return (cos(eta*eps))*eval_S_T(eps,params);
  //return (1.)*eval_S_T(eps,params);
}

double _SL_cos_integrand(double eps, void *params) {
  //double eps = - (1.-1/X);
  double eta = ((double *)params)[0];
  //double eps = t/eta;
  //cout << " eta = " << eta << " , eps = " << eps << endl;
  return (cos(eta*eps))*eval_S_L(eps,params);
}

/*--------------------------------------------------------------------*/
// define the special functions g and h

double h_T(double eta) {
  double params[1];
  double a=0.;
  double res, err;
  params[0] = eta;
  //cout << " hT:" << endl;
  integrate_osc(a,eval_S_T,params,&res,&err,GSL_INTEG_SINE);
  return 2.*res;
}

double h_L(double eta) {
  double params[1];
  double a=0.;
  double res, err;
  params[0] = eta;
  //cout << " hL:" << endl;
  integrate_osc(a,eval_S_L,params,&res,&err,GSL_INTEG_SINE);
  return 2.*res;
}

double g_T(double eta) {
  double params[1];
  params[0] = eta;
  double a=M_PI/2./eta;//0.;
  double res, err;
  double resT2, res3;
  double eps_max = 10.;

  params[0] = eta;
  double temp=0.;

  //cout << " gT: [a,inf]" << endl;
  integrator(0.,a,_ST_cos_integrand,params,&res3,&err);
  integrate_osc(a,eval_S_T,params,&resT2,&err,GSL_INTEG_COSINE);
  return 2.*(cT-res3-resT2);
}

double g_L(double eta) {
  double params[1];
  double a=M_PI/2./eta;//0.;
  double res, err;
  double resL2, res3;
  params[0] = eta;
  //cout << " gL: " << endl;
  integrate_osc(0.,eval_S_L,params,&resL2,&err,GSL_INTEG_COSINE);
  return 2.*(cL-resL2);
}

/*--------------------------------------------------------------------*/
// the function chi(eta,mD/T)

double nB(double x) {
  if (x>0) {
    double e = exp(-x), denom;
    denom = - expm1(-x);
    return e/denom;
  }
  else { return -(1.+nB(-x)); }
}

double _chi_integrand(double eps, void *params) {
  double eta       = ((double *)params)[0];
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return (1.-cos(eta*eps))*res;
}

double _chi_one_integrand(double X, void *params) {
  double eta       = ((double *)params)[0];
  double mD_over_T = ((double *)params)[1];
  double a         = ((double *)params)[2];
  double eps = a-(1.-1./X);
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res/sqr(X);
}

double _chi_cos_integrand(double eps, void *params) {
  double eta       = ((double *)params)[0];
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res; // omit cos(eta*eps) factor!
}

double chi(double eta, double mD_over_T) {
  if (mD_over_T < 1e-3) { return 0.; }
  double a=M_PI/2./eta;//0.;
  double params[3];
  double res1, res2, res3, err;
  params[0] = eta;
  params[1] = mD_over_T;
  params[2] = a;
  //cout << " gL: " << endl;
  integrator(0.,a,_chi_integrand,params,&res1,&err);
  integrator(0.,1.,_chi_one_integrand,params,&res2,&err);
  integrate_osc(a,_chi_cos_integrand,params,&res3,&err,GSL_INTEG_COSINE);
  return 4.*(res1+res2-res3);
}

double _chi_A(double eps, void *params) {
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res - 1./(3.*eps*mD_over_T);
}

double _chi_B(double X, void *params) {
  double eps = 1./X;
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res/sqr(X);
}

double _chi_D(double X, void *params) {
  double eps = - 1. + 1./X;
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res*sqr(eps)/sqr(X);

}

double chi_large_eta_const(double mD_over_T) {
  double params[2];
  double res1, res2, err;
  params[0] = 0.;
  params[1] = mD_over_T;
  integrator(0.,1.,_chi_A,params,&res1,&err);
  integrator(0.,1.,_chi_B,params,&res2,&err);
  return res1+res2;
}

double chi_small_eta_const(double mD_over_T) {
  double params[2];
  double res, err;
  params[0] = 0.;
  params[1] = mD_over_T;
  integrator(1e-4,1.,_chi_D,params,&res,&err);
  return 2.*res;
}

/*--------------------------------------------------------------------*/

int main() {

  cout.precision(8);
    double eta_=0.01, mD_over_T=.3;
  double res = chi(eta_,mD_over_T);
  double cst = chi_large_eta_const(mD_over_T);
  double c2  = chi_small_eta_const(mD_over_T);
  cout << " eta  = " << eta_ << endl;
  cout << " mD/T = " << mD_over_T << endl;
  cout << " chi  = " << res << endl;
  cout << " cst  = " << cst << endl;
  cout << " est. = " << 4.*( (1./(3.*mD_over_T))*( log(eta_) + GAMMA_E ) + cst  ) << endl;
  cout << " c2= " << c2<< endl;
  cout << " small= " << c2*sqr(eta_) << endl;
  return 0;

  double s, s_err;
  double params[1] = {0.};

  cout << "\n --> checking normalisation: " << endl;

  integrator(0.,1.,_S_integrand,params,&s,&s_err);
  cout << " c0  = " << s << " [err = " << s_err << "]"<< endl;

  integrator(0.,1.,_ST_integrand,params,&s,&s_err);
  cout << " cT  = " << s << " [err = " << s_err << "]"<< endl;

  integrator(0.,1.,_SL_integrand,params,&s,&s_err);
  cout << " cL  = " << s << " [err = " << s_err << "]" << endl;

  cout << "\n --> checking constants in large-x limit: " << endl;

  double temp;
  integrator(0.,1.,_S_times_eps,params,&s,&s_err);
  cout << " c1  = " << s << " [err = " << s_err << "]"<< endl;
  temp += s;

  integrator(0.,1.,_S_subtr,params,&s,&s_err);
  cout << " c2  = " << s << " [err = " << s_err << "]"<< endl;
  temp += s;

  temp *= 2.;
  temp += 1.-GAMMA_E;
  cout << " 2(c1+c2)+1-gamma  = " << temp << endl;

  return 0;
  /*

  cout << "\n --> tabulating g & h functions: " << endl;

  fout.open("data/table_g_h.dat");
  fout.precision(8);
  fout << "# columns:\n";
  fout << "# eta,            hT,               hL,               gT,               gL\n";

  double eta, hT, hL, gT, gL;
  double eta_min, eta_max;

  int eta_N = 5000;

  eta_min = 1e-3, eta_max=1e1;
  double delta  = (eta_max-eta_min)/((double)eta_N-1);
  eta = eta_min;

  loop(i,0,eta_N) { // first loop is linear scale from eta ~ [0,10]

    if (eta>1e2) { tolosc = 1e-7; }
    hT = h_T(eta);
    hL = h_L(eta);
    gT = g_T(eta);
    gL = g_L(eta);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << hT
         <<     "    " << hL
         <<     "    " << gT
         <<     "    " << gL
         << endl;

    eta += delta;
  }

  eta_min = 1e1, eta_max=1e3;
  double ratio  = pow(eta_max/eta_min,1./((double)eta_N-1));
  eta = eta_min;

  loop(i,0,eta_N) { // second loop is log scale for eta ~ [10,1000]

    if (eta>1e2) { tolosc = 1e-7; }
    hT = h_T(eta);
    hL = h_L(eta);
    gT = g_T(eta);
    gL = g_L(eta);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << hT
         <<     "    " << hL
         <<     "    " << gT
         <<     "    " << gL
         << endl;

    eta *= ratio;
  }

  fout.close();
  cout << " finished... " << endl;//*/

  cout << "\n --> tabulating the chi function: " << endl;

  fout.open("data/table_chi.dat");
  fout.precision(8);
  fout << "# columns:\n";
  fout << "# eta,            chi(0.1),         chi(0.3),         chi(1.0),         chi(10.)\n";

  double eta, chi_a, chi_b, chi_c, chi_d, chi_e;
  double eta_min, eta_max;

  int eta_N = 5000;

  eta_min = 1e-3, eta_max=1e1;
  double delta  = (eta_max-eta_min)/((double)eta_N-1);
  eta = eta_min;

  loop(i,0,eta_N) { // first loop is linear scale from eta ~ [0,10]

    if (eta>1e2) { tolosc = 1e-6; }
    chi_a = chi(eta,.1);
    chi_b = chi(eta,.3);
    chi_c = chi(eta,1.);
    chi_d = chi(eta,10.);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << chi_a
         <<     "    " << chi_b
         <<     "    " << chi_c
         <<     "    " << chi_d
         << endl;

    eta += delta;
  }

  eta_min = 1e1, eta_max=1e3;
  double ratio  = pow(eta_max/eta_min,1./((double)eta_N-1));
  eta = eta_min;

  loop(i,0,eta_N) { // second loop is log scale for eta ~ [10,1000]

    if (eta>1e2) { tolosc = 1e-6; }
    chi_a = chi(eta,.1);
    chi_b = chi(eta,.3);
    chi_c = chi(eta,1.);
    chi_d = chi(eta,10.);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << chi_a
         <<     "    " << chi_b
         <<     "    " << chi_c
         <<     "    " << chi_d
         << endl;

    eta *= ratio;
  }

  fout.close();
  cout << " finished... " << endl;//*/


  return 0;
}
