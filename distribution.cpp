#include "inc/macros.hpp"
#include "inc/integrators.hpp"
#include "inc/htl.hpp"
#include <gsl/gsl_sf_gamma.h>

using namespace std;

double const c1  =  0.143575;
double const c2  = -0.123714;

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
  if (T<.01) { 
    //cout << " [T=0] case triggered ! " << endl;
    if (eps<mu) { return 1.-sqr(eps/mu)/6.; }
    else { return .5+mu/eps/3.; }
  }

  if (eps<1e-3) { return 1. - 3.*T*eps/(2.*sqr(M_PI)*mD2) - (nf-3.)*sqr(eps)/(12.*sqr(M_PI)*mD2); }

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
  double eps = a-(1.-1./X);
  double res = eval_S_full(eps,params)*eps - 1/(4.*eps);
  return res/sqr(X);
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

double _kappa(double T, double mu, double nf) {
  double res, err;
  double params[4] = {0.,T,mu,nf};

  double ka=0.;
  double Lambda_star=1.; // = arbitrary separation scale
  cout << " Lambda_star = " << Lambda_star << endl;
  params[0] = Lambda_star;
  // --> integral of S_full(eps)*eps from eps=[0,Lambda_star]
  integrator(0.,1.,_S_full_times_eps,params,&res,&err);
  cout << " a1  = " << res << " [err = " << err << "]"<< endl;
  ka += res;

  // --> integral eps*[ S_full(eps) - 1/(4*eps^2) ] from [Lambda_star,+inf)
  integrator(0.,1.,_S_full_subtr,params,&res,&err);
  cout << " a2  = " << res << " [err = " << err << "]"<< endl;
  ka += res;

  ka *= 4.;
  ka += 1.-GAMMA_E;
  ka -= log(Lambda_star);
  cout << " 4(a1+a2)+1-gamma-log(Lam/mD)  = " << ka << endl; 
  
  return ka;
}



/*--------------------------------------------------------------------*/
// convergence acceleration

struct Levin {
  vd num, den;
  int n, ncv;
  bool conv; // convergence flag
  double small, big;
  double eps, lastval, lasteps;

  Levin(int nmax, double _eps) : num(nmax), den(nmax), n(0), ncv(0),
                                 conv(false), eps(_eps), lastval(0.) {
    small = numeric_limits<double>::min()*10.;
    big   = numeric_limits<double>::max();
  }

  double next(double sum, double omega, double beta=1.) {
    double fact, ratio, term, val;
    term   = 1./(beta+n);
    den[n] = term/omega;
    num[n] = sum*den[n];
    if (n > 0) {
      ratio = (beta+n-1)*term;
      loop(j,1,n+1) {
        fact = (n-j+beta)*term;
        num[n-j] = num[n-j+1] - fact*num[n-j];
        den[n-j] = den[n-j+1] - fact*den[n-j];
        term *= ratio;
      }
    }
    n++;
    val = fabs(den[0]) < small ? lastval : num[0]/den[0];
    lasteps = fabs(val-lastval);
    if (lasteps <= eps) ncv++;
    if (ncv >= 2) conv = 1;
    return (lastval = val);
  }
};

/*--------------------------------------------------------------------*/
#include <gsl/gsl_spline.h>
// for interpolations:

namespace itp {

  static gsl_interp_accel *accelerator = nullptr;
  static gsl_spline    *spline_H = nullptr;
  static gsl_spline    *spline_G = nullptr;
  static gsl_spline *spline_DelG = nullptr;

  double T, mu, nf;
  double r_large;
  double r_small;
  double a, b, d;
  double c = (16./9.)*pow(2./sqr(M_PI),2./3.)*cos(5*M_PI/6.)*gsl_sf_gamma(5./3.);

  void init(double _T, double _mu, double _nf) {

    cout << dots << bold << "Preparing interpolation table..." << reset << endl;
    vd eta_list, G_list, DelG_list, H_list; // lists for interpolation

    // set main parameters:
    T  = _T;
    mu = _mu;
    nf = _nf;
    cout << "parameters (in units of mD) --> T  = " << T  << " (temperature)\n"
         << "                                mu = " << mu << " (chemical potential)\n"
         << "                                nf = " << nf << " (number of light flavours)\n";

    // some constants needed later:
    a = DeltaG_small_eta_const(T,mu,nf);
    b = G_large_eta_const(T,mu,nf);
    d = H_small_eta_const(T,mu,nf);

    // filenaming convention:
    stringstream filename;
    filename << "data/table_G_H_{T=" << fixed << setprecision(2) << T << ",mu=" << mu << ",nf=" << nf << "}.dat";

    fin.open(filename.str());

    // TODO: allow option to regenerate data if file doesn't exist!
    if (!fin.is_open()) {
      cerr << "Error: [" << filename.str() << "] must be generated first!\n";
      exit(EXIT_FAILURE);
    }

    cout << "input file (" << filename.str() << ") successfully loaded...\n";

    string line;
    while (getline(fin, line)) {
      if (line.empty() || line[0]=='#') continue;
      stringstream ss(line);
      double _eta, _G, _DelG, _H;
      ss >> _eta >> _G >> _DelG >> _H;
      eta_list.push_back(_eta);
      G_list.push_back(_G);
      DelG_list.push_back(_DelG);
      H_list.push_back(_H);
    }

    // generic constructor for GSL interpolation:
    auto mk_spline = [](vd &x, vd &y) {
      gsl_spline *spl = gsl_spline_alloc(gsl_interp_cspline, x.sz);
      gsl_spline_init(spl,x.ar,y.ar,x.sz);
      return spl;
    };

    accelerator  = gsl_interp_accel_alloc();
    spline_G     = mk_spline(eta_list,G_list);
    spline_DelG  = mk_spline(eta_list,DelG_list);
    spline_H     = mk_spline(eta_list,H_list);

    cout << "--> interpolation complete!\n";
    fin.close();
  }

  void free() {
    gsl_spline_free(spline_G);
    gsl_spline_free(spline_DelG);
    gsl_spline_free(spline_H);
    gsl_interp_accel_free(accelerator);
  }

  double G_func(double eta) {
    if (eta<1e-3) { return eta*M_PI/4.; }
    else if (eta>1e3) { return b + c*pow(eta,-5./3.); }
    else { return gsl_spline_eval(spline_G, eta, accelerator); }
  }
  double DelG_func(double eta) {
    if (T<1e-2) { return 0.; }
    else if (eta<1e-3) { return a*sqr(eta); }
    else if (eta>1e3) { return 4.*T/3.*( log(eta) + GAMMA_E ); }
    else { return gsl_spline_eval(spline_DelG, eta, accelerator); }
  }
  double H_func(double eta) {
    if (eta<1e-3) { return .5*eta*( log(1./eta) + 4.*d - GAMMA_E + 1. ); }
    else if (eta>1e3) { return 2./(3.*eta); }
    else { return gsl_spline_eval(spline_H, eta, accelerator); }
  }
  
  double _eta_integrand(double eta, void *params) {
    double x = ((double *)params)[0];
    double Del = ((double *)params)[1];
    double G = G_func(eta)+DelG_func(eta);
    double H = H_func(eta);
    double b = G_large_eta_const(T,mu,nf);
    //double res = 
    //cos( eta*Del - x*h )*exp( -x*(g-2.*c0) ) - cos( eta*Del );
    //return res*exp( -x*2.*c0 ) ;
    return cos( eta*Del - x*2.*H )*exp( -x*2.*G ) - exp(-x*2.*b)*cos( eta*Del );
  }

  double lev_int(double x, double Del) {
    double ans_0, err_0;
    double params[2] = {x,Del};
  //if (Del<1.) {
  //  integrator(0.01,1.,_eta_integrand2,params,&ans_0,&err_0);
  //  return ans_0/M_PI;
  //} else {
    double width = M_PI/max( fabs(Del), 1e-3 ); // safety padding if Del=0
    int nterm = 100;
    double beta=.5, a=width/2., b=width/2., sum=.0;
    Levin series(100,0.);
    cout << setw(5) << "N" << setw(19) << "Sum (direct)" << setw(21)
         << "Sum (Levin)" << endl;
    integrator(0.,width/2.,_eta_integrand,params,&ans_0,&err_0);
    sum += ans_0;

    double ans, prev;
    loop(n,0,nterm+1) {
      prev = ans;
      b += width;
      double s, s_err;
      integrator(a,b,_eta_integrand,params,&s,&s_err);
      a = b;
      sum += s;
      //if (fabs(s)<epsabs) { continue; }
      double omega = (beta+n)*s;
      ans   = series.next(sum,omega,beta);
      if (isnan(ans)) { 
        cout << " Del = " << Del << " , ans_0 = " << ans_0 << " , width = " << width << endl; 
        cout << " a = " << a << " , b = " << b << " , s = " << s << endl << endl; 
      }
      cout << setw(5) << n << fixed << setprecision(14) << setw(21)
          << sum << setw(21) << ans << endl;
      if ((fabs(ans-prev)<epsabs)||(n>50)) { return ans/M_PI; }
    }
    return ans/M_PI;
  }

  // file writing (usual parameters)
  void scan_Del(double xbar, string fname) {

    fout.open(fname);
    fout.precision(8);

    fout << "# columns: Delta/mD, mD*f(xbar,Delta)" << endl;

    double f;
    double Del = .01;
    Del = -10.05013;
    Del = -1.;


    //while (Del < 30.) {
    while (Del < 3.) {
      f = lev_int(xbar,Del);
      fout << scientific << Del << "    " << f << endl;
      Del += .05;
    }

    fout.close();
    cout << " finished... " << endl;//*/
                                  //
  }

  // file writing (for small Delta)
  void scan_small_D(double xbar, string fname) {

    fout.open(fname);
    fout.precision(8);

    fout << "# columns: Delta/T, mD*f(Del>0), mD*f(Del<0), approx." << endl;

    double f_pos, f_neg, f_app;
    double Del = .01; // = Delta/T
    //Del = -1.5013;


    while (Del < 10.) {
      f_pos = lev_int(xbar,+Del*T);
      f_neg = lev_int(xbar,-Del*T);
      f_app = (2./3.)*xbar*exp(-2.*c0*xbar)*pow( fabs(Del), -1.+(4./3.)*xbar*T );
      fout << scientific << Del << "    " << f_pos
                                << "    " << f_neg
                                << "    " << f_app << endl;
      //Del += .01;
      Del *= 1.01;
    }

    fout.close();
    cout << " finished... " << endl;//*/
                                  //
  }


  // file writing (i.t.o. rescaled variables)
  void scan_scaled(double xbar, string fname) {

    fout.open(fname);
    fout.precision(8);
    //double C = 2*(c1+c2) + 1. - GAMMA_E;
    cout << dots << bold << "Starting scan for {T=" << T << ",mu=" << mu << ",nf=" << nf << "}" << reset << endl;
    double kapp =  _kappa(T,mu,nf);
    cout << "    --> kappa = " << kapp << endl;

    fout << "# columns: Delta/mD/xbar-log(xbar)-C, xbar*mD*f(xbar,Delta)" << endl;

    double f;
    double Del_over_xbar = .01;
    Del_over_xbar = -5.0 + log(xbar) + kapp;


    while (Del_over_xbar-log(xbar)-kapp<15) {
      f = lev_int(xbar,Del_over_xbar*xbar)*xbar;
      fout << scientific << Del_over_xbar-log(xbar)-kapp << "    " << f << endl;
      Del_over_xbar += .05;
    }

    fout.close();
    cout << " finished... " << endl;//*/
                                  //
  }


}

/*--------------------------------------------------------------------*/

bool askToProceed() {
  string input;
  cout << dots << bold << "Proceed with scan? [Y/n] " << reset;
  getline(cin, input);
  if (input.empty()||(input=="Y")||(input == "y")) return true;
  if ((input == "N")||(input == "n")) return false;
  return false;
}

/*--------------------------------------------------------------------*/

int main() {
  using namespace itp;

  init(.1,10.,3.);

  if (!askToProceed()) return 0;

  scan_Del(.2,"data/f_T_0p1_mu_10_x_0p2.dat");

  //double eta = .00101;
  //cout << " eta   = " << eta          << " :\n";
  //cout << " G     = " << G_func(eta)  << endl;
  //cout << " DelG  = " << DelG_func(eta) << endl;
  //cout << " H     = " << H_func(eta)  << endl;

  //eta = .00099;
  //cout << " eta   = " << eta          << " :\n";
  //cout << " G     = " << G_func(eta)  << endl;
  //cout << " DelG  = " << DelG_func(eta) << endl;
  //cout << " H     = " << H_func(eta)  << endl;



  /*
  cout << " R(eps=0.1,T=1,nf=3) = " << R_hard(.1,1.,0.,3.) << endl;
  cout << " R(eps=10,T=1,nf=3) = " << R_hard(10.,1.,0.,3.) << endl;
  cout << " R(eps=100,T=1,nf=3) = " << R_hard(100.,1.,0.,3.) << endl;

  R_fixed_T_mu(1.0,0,"data/calR_T1_mu0.dat"); 
  R_fixed_T_mu(10.,0,"data/calR_T10_mu0.dat"); 
  R_fixed_T_mu(0.01,1.0,"data/calR_T0_mu1.dat"); 
  R_fixed_T_mu(0.01,10.,"data/calR_T0_mu10.dat"); 
  R_fixed_T_mu(3.,1.,"data/calR_T3_mu1.dat"); //*/

  //cout << " S(eps=0.1,T=1,nf=3) = " << S_full(.1,1.,0.,3.) << endl;
  //cout << " S(eps=10,T=1,nf=3) = " << S_full(10.,1.,0.,3.) << endl;
  //cout << " S(eps=100,T=1,nf=3) = " << S_full(100.,1.,0.,3.) << endl;

  //scan_S_full(1./3.,0,"data/S_full_T0p3_mu0.dat"); 
  //scan_S_full(1.0,0,"data/S_full_T1_mu0.dat"); 
  //scan_S_full(3.0,0,"data/S_full_T3_mu0.dat"); 
  //scan_S_full(10.0,0,"data/S_full_T10_mu0.dat"); 

/*
  scan_Del(.2,"data/f_T_0p33_mu_0_x_0p2.dat");
  scan_Del(.4,"data/f_T_0p33_mu_0_x_0p4.dat");
  scan_Del(.6,"data/f_T_0p33_mu_0_x_0p6.dat");
  scan_Del(.8,"data/f_T_0p33_mu_0_x_0p8.dat");
  scan_Del(1.,"data/f_T_0p33_mu_0_x_1p0.dat");
  scan_Del(2.,"data/f_T_0p33_mu_0_x_2p0.dat");
  scan_Del(4.,"data/f_T_0p33_mu_0_x_4p0.dat");
  scan_Del(8.,"data/f_T_0p33_mu_0_x_8p0.dat");
  scan_Del(16.,"data/f_T_0p33_mu_0_x_16p0.dat");//*/



  //double ka = _kappa(T,mu,nf);
 
 /* 
  scan_scaled(.2,"data/phi_T_0p33_mu_0_x_0p2.dat");
  scan_scaled(.4,"data/phi_T_0p33_mu_0_x_0p4.dat");
  scan_scaled(.6,"data/phi_T_0p33_mu_0_x_0p6.dat");
  scan_scaled(.8,"data/phi_T_0p33_mu_0_x_0p8.dat");
  scan_scaled(1.,"data/phi_T_0p33_mu_0_x_1p0.dat");
  scan_scaled(2.,"data/phi_T_0p33_mu_0_x_2p0.dat");
  scan_scaled(4.,"data/phi_T_0p33_mu_0_x_4p0.dat");
  scan_scaled(6.,"data/phi_T_0p33_mu_0_x_6p0.dat");
  scan_scaled(8.,"data/phi_T_0p33_mu_0_x_8p0.dat");
  scan_scaled(16.,"data/phi_T_0p33_mu_0_x_16p0.dat");//*/


  //double eta = .5;
  //
  //cout << " eta = " << eta      << " :\n";
  //cout << " hT  = " << h_T(eta) << endl;
  //cout << " hL  = " << h_L(eta) << endl;
  //cout << " gT  = " << g_T(eta) << endl;
  ////cout << " gL  = " << g_L(eta) << endl;
  //cout << " r   = " <<   r(eta) << endl;

  //scan_small_D(.6,"data/f_zoom_kappa_0p5_x_0p6.dat");
  //scan_Del(.6,"data/f_kappa_0p5_x_0p6.dat");
  /*
  scan_Del(.2,"data/f_kappa_1p0_x_0p2.dat");
  scan_Del(.4,"data/f_kappa_1p0_x_0p4.dat");
  scan_Del(.6,"data/f_kappa_1p0_x_0p6.dat");
  scan_Del(.8,"data/f_kappa_1p0_x_0p8.dat");
  scan_Del(2.,"data/f_kappa_1p0_x_2p0.dat");
  scan_Del(8.,"data/f_kappa_1p0_x_8p0.dat");

  scan_scaled(.2,"data/scaled_f_kappa_1p0_x_0p2.dat");
  scan_scaled(.4,"data/scaled_f_kappa_1p0_x_0p4.dat");
  scan_scaled(.6,"data/scaled_f_kappa_1p0_x_0p6.dat");
  scan_scaled(.8,"data/scaled_f_kappa_1p0_x_0p8.dat");
  scan_scaled(2.,"data/scaled_f_kappa_1p0_x_2p0.dat");
  scan_scaled(8.,"data/scaled_f_kappa_1p0_x_8p0.dat");//*/

  
  //scan_scaled(4.,"data/scaled_f_T_3mD_x_4p0.dat");
  //scan_Del(4.,"data/f_T_3mD_x_4p0.dat");
  /*
  scan_Del(.2,"data/f_T_3mD_x_0p2.dat");
  scan_Del(.4,"data/f_T_3mD_x_0p4.dat");
  scan_Del(.6,"data/f_T_3mD_x_0p6.dat");
  scan_Del(.8,"data/f_T_3mD_x_0p8.dat");
  scan_Del(2.,"data/f_T_3mD_x_2p0.dat");
  scan_Del(8.,"data/f_T_3mD_x_8p0.dat");

  scan_scaled(.2,"data/scaled_f_T_3mD_x_0p2.dat");
  scan_scaled(.4,"data/scaled_f_T_3mD_x_0p4.dat");
  scan_scaled(.6,"data/scaled_f_T_3mD_x_0p6.dat");
  scan_scaled(.8,"data/scaled_f_T_3mD_x_0p8.dat");
  scan_scaled(2.,"data/scaled_f_T_3mD_x_2p0.dat");
  scan_scaled(8.,"data/scaled_f_T_3mD_x_8p0.dat");//*/

/*
  scan_Del(.2,"data/f_T_0p3mD_x_0p2.dat");
  scan_Del(.4,"data/f_T_0p3mD_x_0p4.dat");
  scan_Del(.6,"data/f_T_0p3mD_x_0p6.dat");
  scan_Del(.8,"data/f_T_0p3mD_x_0p8.dat");
  scan_Del(1.,"data/f_T_0p3mD_x_1p0.dat");
  scan_Del(2.,"data/f_T_0p3mD_x_2p0.dat");
  scan_Del(4.,"data/f_T_0p3mD_x_4p0.dat");
  scan_Del(8.,"data/f_T_0p3mD_x_8p0.dat");/*

  scan_scaled(.2,"data/scaled_f_kappa_0p3_x_0p2.dat");
  scan_scaled(.4,"data/scaled_f_kappa_0p3_x_0p4.dat");
  scan_scaled(.6,"data/scaled_f_kappa_0p3_x_0p6.dat");
  scan_scaled(.8,"data/scaled_f_kappa_0p3_x_0p8.dat");
  scan_scaled(2.,"data/scaled_f_kappa_0p3_x_2p0.dat");
  scan_scaled(8.,"data/scaled_f_kappa_0p3_x_8p0.dat");//*/

/*
  scan_Del(.2,"data/f_kappa_10_x_0p2.dat");
  scan_Del(.4,"data/f_kappa_10_x_0p4.dat");
  scan_Del(.6,"data/f_kappa_10_x_0p6.dat");
  scan_Del(.8,"data/f_kappa_10_x_0p8.dat");
  scan_Del(2.,"data/f_kappa_10_x_2p0.dat");
  scan_Del(8.,"data/f_kappa_10_x_8p0.dat");

  scan_scaled(.2,"data/scaled_f_kappa_10_x_0p2.dat");
  scan_scaled(.4,"data/scaled_f_kappa_10_x_0p4.dat");
  scan_scaled(.6,"data/scaled_f_kappa_10_x_0p6.dat");
  scan_scaled(.8,"data/scaled_f_kappa_10_x_0p8.dat");
  scan_scaled(2.,"data/scaled_f_kappa_10_x_2p0.dat");
  scan_scaled(8.,"data/scaled_f_kappa_10_x_8p0.dat");//*/


  return 0;
}
