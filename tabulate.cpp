#include "inc/macros.hpp"
#include "inc/integrators.hpp"
#include "inc/htl.hpp"

using namespace std;

/*--------------------------------------------------------------------*/
// hard part and interpolation

#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>

double Li2(double z) { return gsl_sf_dilog(z); }

double Li3(double z) { // only for z=(-\infty,1]
  double res = 0., temp, Z3 = gsl_sf_zeta(3), lz = log(fabs(z));

  if (z<-1.)        { return Li3(1./z) - lz*( sqr(lz) + sqr(M_PI) )/6. ; }
  else if (z==1.)   { return      Z3; }
  else if (z==-1.)  { return -.75*Z3; }
  else if (z<0.)    { return .25*Li3(sqr(z))-Li3(-z); }
  else if (z<.25)   {
    temp = z;
    for (int i=1;i<10;i++) {
      res  += temp/cube((double)i);
      temp *= z;
    }
  }
  else if (z>.25)   {
    res   = Z3 + sqr(M_PI)*lz/6.+.5*sqr(lz)*(1.5-log(-lz));
    temp  = 2.;
    for (int i=3;i<10;i++) {
      temp *= (double)i;
      res+= gsl_sf_zeta(3-i)*pow(lz,i)/temp;
    }
  }
  return res;
}

double l1f(double x) { return log(1+exp(-x)); }
double l2f(double x) { 
  if (x<-50.) {
    return -.5*sqr(x) - sqr(M_PI)/6.;
  }
  return Li2(-exp(-x)); 
}
double l3f(double x) { 
  if (x<-50.) {
    return sqr(x)*(x+sqr(M_PI))/6.;
  }
  return Li3(-exp(-x)) ; 
}
double l1b(double x) { return log(-expm1(-x)); }
double l2b(double x) { 
  if (x>50.) {
    return exp(-x);
  }
  return Li2(exp(-x)) ; 
}
double l3b(double x) { 
  if (x>50.) {
    return exp(-x);
  }
  return Li3(exp(-x)) ; 
}

double R_hard(double eps, double T, double mu, double nf) {
  double mD2 = sqr(mu)*nf/(2.*sqr(M_PI)) + sqr(T)*(1.+nf/6.); // factor of g^2 not included!

  //if (eps<1e-3) { return 1.; }
  if (eps<1e-3) { return 1. - 3.*T*eps/(2.*sqr(M_PI)*mD2) - (nf-3.)*sqr(eps)/(12.*sqr(M_PI)*mD2); }
  if (T<.1) { 
    //cout << " [T=0] case triggered ! " << endl;
    if (eps<mu) { return 1.-sqr(eps/mu)/6.; }
    else { return .5+mu/eps/3.; }
  }

  double Phi_f, Phi_b;

  Phi_b = eps*sqr(T)*( sqr(M_PI)/6. - l2b(eps/T) ) + 2.*cube(T)*( gsl_sf_zeta(3.)-l3b(eps/T) );

  Phi_f  = eps*sqr(T)*( l2f((eps+mu)/T)-l2f(+mu/T) ) + 2.*cube(T)*( l3f((eps+mu)/T)-l3f(+mu/T) );
  Phi_f += eps*sqr(T)*( l2f((eps-mu)/T)-l2f(-mu/T) ) + 2.*cube(T)*( l3f((eps-mu)/T)-l3f(-mu/T) );

  double res = 6.*Phi_b + nf*Phi_f;
  res /= cube(eps);
  res /= sqr(2.*M_PI);
  res /= mD2;
  res *= 2*sqr(eps);
  return res;
}

// file writing (i.t.o. rescaled variables)
void R_fixed_T_mu(double T, double mu, string fname) {

  fout.open(fname);
  fout.precision(8);

  fout << "# columns: eps/max(pi.T,mu), R(nf=0), R(nf=2), R(nf=3)" << endl;

  double calR_1, calR_2, calR_3;
  double eps = 1e-3;
  double Lambda = max(M_PI*T,mu);

  while (eps<1e3) {
    calR_1 = R_hard(eps*Lambda,T,mu,0);
    calR_2 = R_hard(eps*Lambda,T,mu,2);
    calR_3 = R_hard(eps*Lambda,T,mu,3);
    fout << scientific << eps    << "    " 
                       << calR_1 << "    " 
                       << calR_2 << "    " 
                       << calR_3 << endl;
    eps*=1.01;
  }

  fout.close();
  cout << " finished... " << endl;//*/
                                  //
}

double S_full(double eps, double T, double mu, double nf) {
   return eval_S_tot(eps,NULL)*R_hard(eps,T,mu,nf);
}

// file writing (i.t.o. rescaled variables)
void scan_S_full(double T, double mu, string fname) {

  fout.open(fname);
  fout.precision(8);

  fout << "# columns: eps/mD, S[nf=3]" << endl;

  double S_;
  double eps = 1e-2;

  while (eps<1e4) {
    S_ = S_full(eps,T,mu,3);
    fout << scientific << eps  << "    " 
                       << S_   << endl;
    eps*=1.01;
  }

  fout.close();
  cout << " finished... " << endl;//*/
                                  //
}


double eval_S_full(double eps, void *params) {
  double T  = ((double *)params)[1];
  double mu = ((double *)params)[2];
  double nf = ((double *)params)[3];
  double calR = R_hard(eps,T,mu,nf);
  //cout << " eps = " << eps << endl;
  //cout << " calR = " << calR << endl;
  return eval_S_tot(eps,NULL)*calR;
}

double _S_full_integrand(double X, void *params) {
  double a = ((double *)params)[0];
  double eps = a - (1.-1/X);
  return eval_S_full(eps,params)/sqr(X);
}

double _S_full_times_eps(double X, void *params) {
  double a = ((double *)params)[0];
  double eps = a*(X);
  return a*eps*eval_S_full(eps,params);
}

double _S_full_subtr(double X, void *params) {
  double a = ((double *)params)[0];
  double eps = a-.1*(1.-1./X);
  //double eps = (X);
  return .1*eps*( eval_S_full(eps,params) - .25/sqr(eps) )/sqr(X);
}

double _S_full_cos_integrand(double eps, void *params) {
  //double eps = - (1.-1/X);
  double eta = ((double *)params)[0];
  //double eps = t/eta;
  //cout << " eta = " << eta << " , eps = " << eps << endl;
  //return (1.-cos(eta*eps))*eval_S_T(eps,params);
  return (cos(eta*eps))*eval_S_full(eps,params);
  //return (1.)*eval_S_T(eps,params);
}



double _H(double eta, double T, double mu, double nf) {
  double params[4];
  double a=0.;
  double res, err;
  params[0] = eta;
  params[1] = T;
  params[2] = mu;
  params[3] = nf;
  //cout << " hT:" << endl;
  integrate_osc(a,eval_S_full,params,&res,&err,GSL_INTEG_SINE);
  return 2.*res;
}

double _G(double eta, double T, double mu, double nf) {
  double params[4];
  //double a=0.;
  double a=M_PI/2./eta;//0.;
  double res1, res2, res3, err;
  params[1] = T;
  params[2] = mu;
  params[3] = nf;
  params[0] = 0.;
  integrator(0.,1.,_S_full_integrand,params,&res1,&err); // S(eps) over [0,+inf)
  params[0] = eta;
  integrator(0.,a,_S_full_cos_integrand,params,&res2,&err); // S(eps)*cos(eta*eps) over [0,a]
  integrate_osc(a,eval_S_full,params,&res3,&err,GSL_INTEG_COSINE); // S(eps)*cos(..) over [a,+inf)
  return 2.*(res1-res2-res3);
}

void tabulate_G_and_H(int eta_N, double T, double mu, double nf) {

  stringstream filename;
  filename << "data/table_G_H_{T=" << fixed << setprecision(2) << T << ",mu=" << mu << ",nf=" << nf << "}.dat";

  fout.open(filename.str());

  fout.precision(8);
  fout << "# columns:\n";
  fout << "# eta,            G,                H\n";

  double eta, G_tmp, H_tmp;
  double eta_min, eta_max;

  //int eta_N = 5000;

  eta_min = 1e-3, eta_max=1e1;
  double delta  = (eta_max-eta_min)/((double)eta_N-1);
  eta = eta_min;

  loop(i,0,eta_N) { // first loop is linear scale from eta ~ [0,10]

    if (eta>1e2) { tolosc = 1e-7; }
    G_tmp = _G(eta,T,mu,nf);
    H_tmp = _H(eta,T,mu,nf);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << G_tmp
         <<     "    " << H_tmp
         << endl;

    eta += delta;
  }

  eta_min = 1e1, eta_max=1e3;
  double ratio  = pow(eta_max/eta_min,1./((double)eta_N-1));
  eta = eta_min*ratio;

  loop(i,1,eta_N) { // second loop is log scale for eta ~ [10,1000]

    if (eta>1e2) { tolosc = 1e-7; }
    G_tmp = _G(eta,T,mu,nf);
    H_tmp = _H(eta,T,mu,nf);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << G_tmp
         <<     "    " << H_tmp
         << endl;

    eta *= ratio;
  }

  fout.close();
  cout << " finished... " << endl;//*/

}



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
// the function r(eta,mD/T)

double nB(double x) {
  if (x>0) {
    double e = exp(-x), denom;
    denom = - expm1(-x);
    return e/denom;
  }
  else { return -(1.+nB(-x)); }
}

double _r_integrand(double eps, void *params) {
  double eta       = ((double *)params)[0];
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return (1.-cos(eta*eps))*res;
}

double _r_one_integrand(double X, void *params) {
  double eta       = ((double *)params)[0];
  double mD_over_T = ((double *)params)[1];
  double a         = ((double *)params)[2];
  double eps = a-(1.-1./X);
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res/sqr(X);
}

double _r_cos_integrand(double eps, void *params) {
  double eta       = ((double *)params)[0];
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res; // omit cos(eta*eps) factor!
}

double r(double eta, double mD_over_T) {
  if (mD_over_T < 1e-3) { return 0.; }
  double a=M_PI/2./eta;//0.;
  double params[3];
  double res1, res2, res3, err;
  params[0] = eta;
  params[1] = mD_over_T;
  params[2] = a;
  //cout << " gL: " << endl;
  integrator(0.,a,_r_integrand,params,&res1,&err);
  integrator(0.,1.,_r_one_integrand,params,&res2,&err);
  integrate_osc(a,_r_cos_integrand,params,&res3,&err,GSL_INTEG_COSINE);
  return 4.*(res1+res2-res3);
}

double _r_A(double eps, void *params) {
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res - 1./(3.*eps*mD_over_T);
}

double _r_B(double X, void *params) {
  double eps = 1./X;
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res/sqr(X);
}

double _r_D(double X, void *params) {
  double eps = - 1. + 1./X;
  double mD_over_T = ((double *)params)[1];
  double res = eval_S_tot(eps,params)*nB(eps*mD_over_T);
  return res*sqr(eps)/sqr(X);

}

double r_large_eta_const(double mD_over_T) {
  double params[2];
  double res1, res2, err;
  params[0] = 0.;
  params[1] = mD_over_T;
  integrator(0.,1.,_r_A,params,&res1,&err);
  integrator(0.,1.,_r_B,params,&res2,&err);
  return res1+res2;
}

double r_small_eta_const(double mD_over_T) {
  double params[2];
  double res, err;
  params[0] = 0.;
  params[1] = mD_over_T;
  integrator(1e-4,1.,_r_D,params,&res,&err);
  return 2.*res;
}

/*--------------------------------------------------------------------*/

void tabulate_g_h(int eta_N) {

  fout.open("data/table_g_h.dat");
  fout.precision(8);
  fout << "# columns:\n";
  fout << "# eta,            hT,               hL,               gT,               gL\n";

  double eta, hT, hL, gT, gL;
  double eta_min, eta_max;

  //int eta_N = 5000;

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
  eta = eta_min*ratio;

  loop(i,1,eta_N) { // second loop is log scale for eta ~ [10,1000]

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

}

void tabulate_r(int eta_N, double mD_over_T) {

  stringstream filename;
  filename << "data/table_r_mDoT=" << fixed << setprecision(2) << mD_over_T << ".dat";

  fout.open(filename.str());
  //fout.open("data/table_r.dat");
  fout << "# columns:\n";
  fout << "# eta,            r(" << fixed << setprecision(2) << mD_over_T << ")\n";
  fout.precision(8);

  double eta, r_a;
  double eta_min, eta_max;

  //int eta_N = 5000;

  eta_min = 1e-3, eta_max=1e1;
  double delta  = (eta_max-eta_min)/((double)eta_N-1);
  eta = eta_min;

  loop(i,0,eta_N) { // first loop is linear scale from eta ~ [0,10]

    if (eta>1e2) { tolosc = 1e-6; }
    r_a = r(eta,mD_over_T);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << r_a
         << endl;

    eta += delta;
  }

  eta_min = 1e1, eta_max=1e3;
  double ratio  = pow(eta_max/eta_min,1./((double)eta_N-1));
  eta = eta_min*ratio;

  loop(i,1,eta_N) { // second loop is log scale for eta ~ [10,1000]

    if (eta>1e2) { tolosc = 1e-6; }
    r_a = r(eta,mD_over_T);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << r_a
         << endl;

    eta *= ratio;
  }

  fout.close();
  cout << " finished... " << endl;//*/


}

/*--------------------------------------------------------------------*/

int main() {

  R_fixed_T_mu(1.0,0,"data/calR_T1_mu0.dat"); 
  R_fixed_T_mu(0,1.0,"data/calR_T0_mu1.dat"); 
  R_fixed_T_mu(1.,3.,"data/calR_T1_mu3.dat"); 
  R_fixed_T_mu(1.,10.,"data/calR_T1_mu10.dat"); //*/

  return 0;
  cout.precision(8);
  double eta_=1., T_over_mD=0., mu_over_mD=10., nf=3.;
  double H_test = _H(eta_,T_over_mD,mu_over_mD,nf);
  double G_test = _G(eta_,T_over_mD,mu_over_mD,nf);
  double hT=h_T(eta_);
  double hL=h_L(eta_);
  double gT=g_T(eta_);
  double gL=g_L(eta_);
  cout << "       H(eta)  = " << H_test << endl;
  cout << "       G(eta)  = " << G_test << endl << endl;

  cout << " OLD: hT(eta)  = " << hT << endl;
  cout << " OLD: hL(eta)  = " << hL << endl;
  cout << " OLD: hL+hL    = " << hL+hT << endl;

  cout << " OLD: gT(eta)  = " << gT << endl;
  cout << " OLD: gL(eta)  = " << gL << endl;
  cout << " OLD: gL+gL    = " << gL+gT << endl << endl;


  double s, s_err;
  double params[4] = {0.,T_over_mD,mu_over_mD,nf};

  cout << "\n --> checking normalisation: " << endl;

  integrator(0.,1.,_S_integrand,params,&s,&s_err);
  cout << " OLD: b   = " << s << " [err = " << s_err << "]"<< endl;

  double temp=0.;
  integrator(0.,1.,_S_times_eps,params,&s,&s_err);
  cout << " a1  = " << s << " [err = " << s_err << "]"<< endl;
  temp += s;

  integrator(0.,1.,_S_subtr,params,&s,&s_err);
  cout << " a2  = " << s << " [err = " << s_err << "]"<< endl;
  temp += s;

  temp *= 2.;
  temp += 1.-GAMMA_E;
  cout << " 2(a1+a2)+1-gamma  = " << temp << endl; //*/
                                                   //
                                                   //
  cout << endl << " NEW: " << endl;
  integrator(0.,1.,_S_full_integrand,params,&s,&s_err);
  cout << "      b   = " << s << " [err = " << s_err << "]"<< endl;


  double kappa=0.;
  double Lambda_star=.4;
  cout << " Lambda_star = " << Lambda_star << endl;
  params[0] = Lambda_star;
  integrator(0.,1.,_S_full_times_eps,params,&s,&s_err);
  cout << " a1  = " << s << " [err = " << s_err << "]"<< endl;
  kappa += s;

  integrator(0.,1.,_S_full_subtr,params,&s,&s_err);
  cout << " a2  = " << s << " [err = " << s_err << "]"<< endl;
  kappa += s;

  kappa *= 4.;
  kappa += 1.-GAMMA_E;
  kappa -= log(Lambda_star);
  cout << " 4(a1+a2)+1-gamma-log(Lam/mD)  = " << kappa << endl; //*/
                                                                //
  kappa=0.;
  Lambda_star=4.3;
  cout << " Lambda_star = " << Lambda_star << endl;
  params[0] = Lambda_star;
  integrator(0.,1.,_S_full_times_eps,params,&s,&s_err);
  cout << " a1  = " << s << " [err = " << s_err << "]"<< endl;
  kappa += s;

  integrator(0.,1.,_S_full_subtr,params,&s,&s_err);
  cout << " a2  = " << s << " [err = " << s_err << "]"<< endl;
  kappa += s;

  kappa *= 4.;
  kappa += 1.-GAMMA_E;
  kappa -= log(Lambda_star);
  cout << " 4(a1+a2)+1-gamma-log(Lam/mD)  = " << kappa << endl; //*/


  //tabulate_G_and_H(5000,T_over_mD,mu_over_mD,nf);
  //tabulate_G_and_H(5000,10.,mu_over_mD,nf);
  //tabulate_G_and_H(5000,1.,mu_over_mD,nf);
  //tabulate_G_and_H(5000,.01,10.,nf);
  //tabulate_G_and_H(5000,3.,1.,nf);

  /*
  double res = r(eta_,mD_over_T);
  double cst = r_large_eta_const(mD_over_T);
  double c2  = r_small_eta_const(mD_over_T);
  cout << " eta  = " << eta_ << endl;
  cout << " mD/T = " << mD_over_T << endl;
  cout << "   r  = " << res << endl;
  cout << " cst  = " << cst << endl;
  cout << " est. = " << 4.*( (1./(3.*mD_over_T))*( log(eta_) + GAMMA_E ) + cst  ) << endl;
  cout << " c2= " << c2<< endl;
  cout << " small= " << c2*sqr(eta_) << endl;

  double s, s_err;
  double params[1] = {0.};

  cout << "\n --> checking normalisation: " << endl;

  integrator(0.,1.,_S_integrand,params,&s,&s_err);
  cout << " b   = " << s << " [err = " << s_err << "]"<< endl;

  integrator(0.,1.,_ST_integrand,params,&s,&s_err);
  cout << " bT  = " << s << " [err = " << s_err << "]"<< endl;

  integrator(0.,1.,_SL_integrand,params,&s,&s_err);
  cout << " bL  = " << s << " [err = " << s_err << "]" << endl;

  cout << "\n --> checking constants in large-x limit: " << endl;

  double temp;
  integrator(0.,1.,_S_times_eps,params,&s,&s_err);
  cout << " a1  = " << s << " [err = " << s_err << "]"<< endl;
  temp += s;

  integrator(0.,1.,_S_subtr,params,&s,&s_err);
  cout << " a2  = " << s << " [err = " << s_err << "]"<< endl;
  temp += s;

  temp *= 2.;
  temp += 1.-GAMMA_E;
  cout << " 2(a1+a2)+1-gamma  = " << temp << endl; //*/

  //return 0;

  //cout << "\n --> tabulating g & h functions: " << endl;

  //tabulate_g_h(5000);

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

  //
  //
  //cout << "\n --> tabulating the r function: " << endl;
  //tabulate_r(5000,0.3);
  //
  //

  /*
  fout.open("data/table_r.dat");
  fout.precision(8);
  fout << "# columns:\n";
  fout << "# eta,            r(0.1),         r(0.3),         r(1.0),         r(10.)\n";

  double eta, r_a, r_b, r_c, r_d, r_e;
  double eta_min, eta_max;

  int eta_N = 5000;

  eta_min = 1e-3, eta_max=1e1;
  double delta  = (eta_max-eta_min)/((double)eta_N-1);
  eta = eta_min;

  loop(i,0,eta_N) { // first loop is linear scale from eta ~ [0,10]

    if (eta>1e2) { tolosc = 1e-6; }
    r_a = r(eta,.1);
    r_b = r(eta,.3);
    r_c = r(eta,1.);
    r_d = r(eta,10.);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << r_a
         <<     "    " << r_b
         <<     "    " << r_c
         <<     "    " << r_d
         << endl;

    eta += delta;
  }

  eta_min = 1e1, eta_max=1e3;
  double ratio  = pow(eta_max/eta_min,1./((double)eta_N-1));
  eta = eta_min;

  loop(i,0,eta_N) { // second loop is log scale for eta ~ [10,1000]

    if (eta>1e2) { tolosc = 1e-6; }
    r_a = r(eta,.1);
    r_b = r(eta,.3);
    r_c = r(eta,1.);
    r_d = r(eta,10.);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << r_a
         <<     "    " << r_b
         <<     "    " << r_c
         <<     "    " << r_d
         << endl;

    eta *= ratio;
  }

  fout.close();
  cout << " finished... " << endl;//*/


  return 0;
}
