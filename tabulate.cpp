#include "inc/macros.hpp"
#include "inc/integrators.hpp"
#include "inc/htl.hpp"

using namespace std;

/*--------------------------------------------------------------------*/
// useful thermal functions:

double nB(double x) {
  if (x>0) {
    double e = exp(-x), denom;
    denom = - expm1(-x);
    return e/denom;
  }
  else { return -(1.+nB(-x)); }
}


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

/*--------------------------------------------------------------------*/
// hard part & interpolation

double R_hard(double eps, double T, double mu, double nf) {
  double mD2 = sqr(mu)*nf/(2.*sqr(M_PI)) + sqr(T)*(1.+nf/6.); // factor of g^2 not included!
                                                              // (in order for T, and mu to 
                                                              // be given "in units of mD", 
                                                              // g should be such that the 
                                                              // (RHS of this eq)*g^2 = 1

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
  /*
   *  S_HTL*R_hard : eps = [a,+inf)
   */
  double a = ((double *)params)[0];
  double eps = a - (1.-1/X);
  return eval_S_full(eps,params)/sqr(X);
}

double _S_full_times_eps(double X, void *params) {
  /*
   *  [ S_HTL*R_hard ]*eps : eps = [0,a]
   */
  double a = ((double *)params)[0];
  double eps = a*(X);
  return a*eps*eval_S_full(eps,params);
}

double _S_full_subtr(double X, void *params) {
  /*
   *  [ S_HTL*R_hard - 1/(4*eps^2) ]*eps : eps = [a,+inf)
   */
  double a = ((double *)params)[0];
  double eps = a-.1*(1.-1./X);
  //double eps = (X);
  return .1*eps*( eval_S_full(eps,params) - .25/sqr(eps) )/sqr(X);
}

double _S_full_cos_integrand(double eps, void *params) {
  /*
   *  S_HTL*R_hard*nB*cos : eps linear
   */
  double eta = ((double *)params)[0];
  //double eps = t/eta;
  //cout << " eta = " << eta << " , eps = " << eps << endl;
  //return (1.-cos(eta*eps))*eval_S_T(eps,params);
  return (cos(eta*eps))*eval_S_full(eps,params);
  //return (1.)*eval_S_T(eps,params);
}

double _DeltaG_integrand(double eps, void *params) {
  /*
   *  S_HTL*R_hard*nB*(1-cos) : eps linear
   */
  double eta   = ((double *)params)[0];
  double T     = ((double *)params)[1];
  double res = eval_S_full(eps,params)*nB(eps/T);
  return (1.-cos(eta*eps))*res;
}

double _DeltaG_one_integrand(double X, void *params) {
  /*
   *  S_HTL*R_hard*nB : eps = [a,+inf)
   */
  double a     = ((double *)params)[0];
  double T     = ((double *)params)[1];
  double eps = a-(1.-1./X);
  double res = eval_S_full(eps,params)*nB(eps/T);
  return res/sqr(X);
}

double _DeltaG_cos_integrand(double eps, void *params) {
  /*
   *  S_HTL*R_hard*nB : eps = linear
   */
  double eta   = ((double *)params)[0];
  double T     = ((double *)params)[1];
  double res = eval_S_full(eps,params)*nB(eps/T);
  return res; // omit cos(eta*eps) factor!
}

//
// large-eta expansion for G
//

double G_large_eta_const(double T, double mu, double nf) {
  double params[4];
  double res, err;
  params[0] = 0.;
  params[1] = T;
  params[2] = mu;
  params[3] = nf;
  integrator(0,1.,_S_full_integrand,params,&res,&err);
  return 2.*res;
}

//
// small-eta expansion for G-tilde
//

double _DeltaG_small(double X, void *params) {
  double eps   = - 1. + 1./X;
  double eta   = ((double *)params)[0];
  double T     = ((double *)params)[1];
  double res = eval_S_full(eps,params)*nB(eps/T);
  return res*sqr(eps/X);
}

double DeltaG_small_eta_const(double T, double mu, double nf) {
  double params[4];
  double res, err;
  params[0] = 0.;
  params[1] = T;
  params[2] = mu;
  params[3] = nf;
  integrator(1e-4,1.,_DeltaG_small,params,&res,&err);
  return 2.*res;
}

//
// small-eta expansion for H:
//

double _H_small_subtr(double X, void *params) {
  double eps   = 1./X;
  double eta   = ((double *)params)[0];
  double T     = ((double *)params)[1];
  double res = eval_S_full(eps,params)*eps - 1/(4.*eps);
  return res/sqr(X);
}

double _H_small_Seps(double eps, void *params) {
  double eta   = ((double *)params)[0];
  double T     = ((double *)params)[1];
  double res = eval_S_full(eps,params)*eps;
  return res;
}

double H_small_eta_const(double T, double mu, double nf) {
  double params[4];
  double res1, res2, err;
  params[0] = 0.;
  params[1] = T;
  params[2] = mu;
  params[3] = nf;
  integrator(0.,1.,_H_small_subtr,params,&res1,&err);
  integrator(0.,1.,_H_small_Seps,params,&res2,&err);
  return res1+res2;
}

/*--------------------------------------------------------------------*/

#include <gsl/gsl_sf_gamma.h>

void print_asympt(double T, double mu, double nf) {
  double a = DeltaG_small_eta_const(T,mu,nf);
  double b = .5*G_large_eta_const(T,mu,nf);
  double small_c = (16./9.)*pow(2./sqr(M_PI),2./3.)*cos(5*M_PI/6.)*gsl_sf_gamma(5./3.);
  double d = H_small_eta_const(T,mu,nf);
  cout << " small-eta expansion: "  << endl;
  cout << "   G(eta) ~ pi*eta/4"    << endl;
  cout << "  DG(eta) ~ a*eta^2" << endl;
  cout << "        a ≈ " << a << endl;
  cout << "   H(eta) ~ eta/2*[ log(1/eta) + 4*d - Gamma_E + 1 ]" << endl;
  cout << "        d ≈ " << d << endl << endl;
  cout << " large-eta expansion: "  << endl;
  cout << "   G(eta) ~ 2*b + c*eta^{-5/3}" << endl;
  cout << "        b ≈ " << b << endl;
  cout << "        c = (16/9)*(2/pi^2)^{2/3}*cos(5pi/6)*Gamma(5/3) ≈ " << small_c << endl;
  cout << "  DG(eta) ~ 4/3*T/mD*[ log(eta) + Gamma_E ]" << endl;
  cout << "   H(eta) ~ 2/(3*eta)" << endl;
}

/*--------------------------------------------------------------------*/

double _DeltaG(double eta, double T, double mu, double nf) {
  if (T < 1e-3) { return 0.; }
  double a=M_PI/2./eta;//0.;
  double params[4];
  double res1, res2, res3, err;
  params[1] = T;
  params[2] = mu;
  params[3] = nf;
  //cout << " gL: " << endl;
  params[0] = eta;
  integrator(0.,a,_DeltaG_integrand,params,&res1,&err);
  integrate_osc(a,_DeltaG_cos_integrand,params,&res3,&err,GSL_INTEG_COSINE);
  params[0] = a;
  integrator(0.,1.,_DeltaG_one_integrand,params,&res2,&err);
  return 4.*(res1+res2-res3);
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
  fout << "# eta,            G,                DeltaG,           H\n";

  double eta, G_tmp, DeltaG_tmp, H_tmp;
  double eta_min, eta_max;

  //int eta_N = 5000;

  eta_min = 1e-3, eta_max=1e1;
  double delta  = (eta_max-eta_min)/((double)eta_N-1);
  eta = eta_min;
  cout << "\n --> beginning tabulation (in eta): N = " << eta_N << " , delta_eta = " << delta << endl;
  cout << "                                    T/mD = " << T << " , mu/mD = " << mu << " , nf = " << nf << endl << endl;

  cout << left
       << setw(12) << " eta: "
       << setw(12) << " G: "
       << setw(12) << " Delta G: "
       << setw(12) << " H: " << endl;

  cout << right << fixed << setprecision(5);

  loop(i,0,eta_N) { // first loop is linear scale from eta ~ [0,10]

    if (eta>1e2) { tolosc = 1e-7; }
    G_tmp  = _G(eta,T,mu,nf);
    DeltaG_tmp = _DeltaG(eta,T,mu,nf);
    H_tmp  = _H(eta,T,mu,nf);
    cout << setw(12) << eta 
         << setw(12) << G_tmp 
         << setw(12) << DeltaG_tmp
         << setw(12) << H_tmp       << endl;

    fout << scientific << eta
         <<     "    " << G_tmp
         <<     "    " << DeltaG_tmp
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
    DeltaG_tmp = _DeltaG(eta,T,mu,nf);
    H_tmp = _H(eta,T,mu,nf);
    cout << eta << endl;

    fout << scientific << eta
         <<     "    " << G_tmp
         <<     "    " << DeltaG_tmp
         <<     "    " << H_tmp
         << endl;

    eta *= ratio;
  }

  fout.close();
  cout << " finished... " << endl;//*/

}

/*--------------------------------------------------------------------*/

double _kappa(double T, double mu, double nf) {
  double res, err;
  double params[4] = {0.,T,mu,nf};

  double kappa=0.;
  double Lambda_star=2.; // = arbitrary separation scale
  //cout << " Lambda_star = " << Lambda_star << endl;
  params[0] = Lambda_star;
  integrator(0.,1.,_S_full_times_eps,params,&res,&err);
  //cout << " a1  = " << res << " [err = " << err << "]"<< endl;
  kappa += res;

  integrator(0.,1.,_S_full_subtr,params,&res,&err);
  //cout << " a2  = " << res << " [err = " << err << "]"<< endl;
  kappa += res;

  kappa *= 4.;
  kappa += 1.-GAMMA_E;
  kappa -= log(Lambda_star);
  //cout << " 4(a1+a2)+1-gamma-log(Lam/mD)  = " << kappa << endl; 
  
  return kappa;
}

/*--------------------------------------------------------------------*/

int main() {

  /*
  R_fixed_T_mu(1.0,0,"data/calR_T1_mu0.dat"); 
  R_fixed_T_mu(0,1.0,"data/calR_T0_mu1.dat"); 
  R_fixed_T_mu(1.,3.,"data/calR_T1_mu3.dat"); 
  R_fixed_T_mu(1.,10.,"data/calR_T1_mu10.dat"); //*/


  cout.precision(8);
  double eta_=1., T_over_mD=20., mu_over_mD=0., nf=3.;
  /*
  double H_test  = _H(eta_,T_over_mD,mu_over_mD,nf);
  double G_test  = _G(eta_,T_over_mD,mu_over_mD,nf);
  double Gt_test = _Gt(eta_,T_over_mD,mu_over_mD,nf);
  double hT=h_T(eta_);
  double hL=h_L(eta_);
  double gT=g_T(eta_);
  double gL=g_L(eta_);
  cout << "       H(eta)  = " << H_test  << endl;
  cout << "       G(eta)  = " << G_test  << endl << endl;
  cout << "      Gt(eta)  = " << Gt_test << endl << endl;

  cout << " OLD: hT(eta)  = " << hT << endl;
  cout << " OLD: hL(eta)  = " << hL << endl;
  cout << " OLD: hL+hL    = " << hL+hT << endl;

  cout << " OLD: gT(eta)  = " << gT << endl;
  cout << " OLD: gL(eta)  = " << gL << endl;
  cout << " OLD: gL+gL    = " << gL+gT << endl << endl;//*/


  double s, s_err;
  double params[4] = {0.,T_over_mD,mu_over_mD,nf};

  cout << "\n --> checking normalisation: " << endl;

  /*
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


  double kapp = _kappa(10.,0.,3.);
  cout << "    kappa ≡ 4(a1+a2)+1-gamma-log(Lam/mD) ≈ " << kapp << endl << endl; 
  /*
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
  cout << " 4(a1+a2)+1-gamma-log(Lam/mD)  = " << kappa << endl; 

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


  print_asympt(T_over_mD,mu_over_mD,nf);
  //double A = Gt_small_eta_const(T_over_mD,mu_over_mD,nf);
  //tabulate_G_and_H(5000,T_over_mD,mu_over_mD,nf);
  //tabulate_G_and_H(5000,10.,mu_over_mD,nf);
  //tabulate_G_and_H(5000,1.,mu_over_mD,nf);
  //tabulate_G_and_H(5000,.0,10.,nf);
  //tabulate_G_and_H(5000,3.,1.,nf);
  //tabulate_G_and_H(5000,1./3.,0.,nf);
  tabulate_G_and_H(5000,2.,10.,nf);

  return 0;
}
