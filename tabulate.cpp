#include "inc/macros.hpp"
#include "inc/integrators.hpp"
#include "inc/htl.hpp"

using namespace std;

const double c0  = 0.72035400;
const double cT  = 0.24944610;
const double cL  = 0.47090790;
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

int main() {

  double s, s_err;
  double params[1] = {0.};

  cout << "\n --> checking normalisation: " << endl;

  integrator(0.,1.,_S_integrand,params,&s,&s_err);
  cout << " c0  = " << s << " [err = " << s_err << "]"<< endl;

  integrator(0.,1.,_ST_integrand,params,&s,&s_err);
  cout << " cT  = " << s << " [err = " << s_err << "]"<< endl;

  integrator(0.,1.,_SL_integrand,params,&s,&s_err);
  cout << " cL  = " << s << " [err = " << s_err << "]" << endl;
  cout << endl;

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

  return 0;
}
