#include "inc/macros.hpp"
#include "inc/integrators.hpp"
#include "inc/htl.hpp"

using namespace std;

double const c1  =  0.143575;
double const c2  = -0.123714;

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
  static gsl_spline *spline_h_T = nullptr;
  static gsl_spline *spline_h_L = nullptr;
  static gsl_spline *spline_g_T = nullptr;
  static gsl_spline *spline_g_L = nullptr;
  static gsl_spline *spline_r   = nullptr;

  double mD_over_T;
  double r_large;
  double r_small;

  void init(double _mD) {
    vd eta_list, hT_list, hL_list, gT_list, gL_list, r_list;
    string line;
    mD_over_T = _mD;

    auto mk_spline = [](vd &x, vd &y) {
      gsl_spline *spl = gsl_spline_alloc(gsl_interp_cspline, x.sz);
      gsl_spline_init(spl,x.ar,y.ar,x.sz);
      return spl;
    };

    //ifstream file("data/table_g_h.dat");
    fin.open("data/table_g_h.dat");

    if (!fin.is_open()) {
      cerr << "Error: [data/table_g_h.dat] must be generated first!\n";
      exit(EXIT_FAILURE);
    }

    while (getline(fin, line)) {
      if (line.empty() || line[0]=='#') continue;
      stringstream ss(line);
      double _eta, _hT, _hL, _gT, _gL;
      ss >> _eta >> _hT >> _hL >> _gT >> _gL;
      eta_list.push_back(_eta);
      hT_list.push_back(_hT);
      hL_list.push_back(_hL);
      gT_list.push_back(_gT);
      gL_list.push_back(_gL);
    }

    accelerator = gsl_interp_accel_alloc();
    spline_h_T  = mk_spline(eta_list,hT_list);
    spline_h_L  = mk_spline(eta_list,hL_list);
    spline_g_T  = mk_spline(eta_list,gT_list);
    spline_g_L  = mk_spline(eta_list,gL_list);

    fin.close();

    if (mD_over_T>0) {
    stringstream filename;
    filename << "data/table_r_mDoT=" << fixed << setprecision(2) << mD_over_T << ".dat";

    // testing limits:
    //r_large = r_large_eta_const(mD_over_T);
    //r_small = r_small_eta_const(mD_over_T);
    //cout << "\n" << " r-large = " << r_large << endl;
    //cout         << " r-small = " << r_small << endl;

    fin.open(filename.str());

    if (!fin.is_open()) {
      cerr << "Error: [" << filename.str() << "] must be generated first!\n";
      exit(EXIT_FAILURE);
    }
    //fin.open("data/table_chi.dat");

    eta_list.clear();
    while (getline(fin, line)) {
      if (line.empty() || line[0]=='#') continue;
      stringstream ss(line);
      double _eta, _r1;// _chi2, _chi3, _chi4;
      ss >> _eta >> _r1;// >> _chi2 >> _chi3 >> _chi4;
      eta_list.push_back(_eta); //(already done)
      r_list.push_back(_r1);
      cout << " eta = " << _eta << " , r = " << _r1 << endl;
    }

    fin.close();

    spline_r  = mk_spline(eta_list,r_list);
    }
  }

  void free() {
    gsl_spline_free(spline_h_T);
    gsl_spline_free(spline_h_L);
    gsl_spline_free(spline_g_T);
    gsl_spline_free(spline_g_L);
    gsl_spline_free(spline_r);
    gsl_interp_accel_free(accelerator);
  }

  double h_T(double eta) {
    if (eta<1e-3) { return (1.-GAMMA_E-log(eta))*eta/3.; }
    else if (eta>1e3) { return 2./(3.*eta); }
    else { return gsl_spline_eval(spline_h_T, eta, accelerator); }
  }
  double h_L(double eta) {
    if (eta<1e-3) { return (1.-GAMMA_E-log(eta))*eta*2./3.; }
    else if (eta>1e3) { return 2.*3.0652/pow(eta,3); }
    else { return gsl_spline_eval(spline_h_L, eta, accelerator); }
  }
  double g_T(double eta) {
    if (eta<1e-3) { return eta*M_PI/6.; }
    else if (eta>1e3) { return 2.*( cT - .239753/pow(eta,5./3.) ); }
    else { return gsl_spline_eval(spline_g_T, eta, accelerator); }
  }
  double g_L(double eta) {
    if (eta<1e-3) { return eta*M_PI/3.; }
    else if (eta>1e3) { return 2.*( cL + .785398/pow(eta,2) ); }
    else { return gsl_spline_eval(spline_g_L, eta, accelerator); }
  }
  double r(double eta) {
    if (mD_over_T<=0.) { return 0.; }
    //if (eta<1e-3) { return (4.7866666)*sqr(eta); } // mD/T = 0.3
    //if (eta<1e-3) { return (0.0014929826)*sqr(eta); } // mD/T = 10.
    if (eta<1e-3) { return (r_small)*sqr(eta); } // mD/T = 10.
    else if (eta>1e3) { 
      //double mD_over_T = .3;
      //double mD_over_T = 10.;
      //return 4.*( (1./(3.*mD_over_T))*( log(eta) + GAMMA_E ) + (-0.046050746) ); // mD/T = 0.3
      //return 4.*( (1./(3.*mD_over_T))*( log(eta) + GAMMA_E ) + (-0.081409316) ); // mD/T = 10.
      return 4.*( (1./(3.*mD_over_T))*( log(eta) + GAMMA_E ) + (r_large) );
    }
    else { return gsl_spline_eval(spline_r, eta, accelerator); }
  }


  double _eta_integrand(double eta, void *params) {
    double x = ((double *)params)[0];
    double Del = ((double *)params)[1];
    double h = h_T(eta) + h_L(eta);
    double g = g_T(eta) + g_L(eta);
    double ch = r(eta);
    //double res = 
    //cos( eta*Del - x*h )*exp( -x*(g-2.*c0) ) - cos( eta*Del );
    //return res*exp( -x*2.*c0 ) ;
    return cos( eta*Del - x*h )*exp( -x*(g+ch) );// - exp(-x*2.*c0)*cos( eta*Del );
  }

  double lev_int(double x, double Del) {
    double ans_0, err_0;
    double params[2] = {x,Del};
  //if (Del<1.) {
  //  integrator(0.01,1.,_eta_integrand2,params,&ans_0,&err_0);
  //  return ans_0/M_PI;
  //} else {
    double width = M_PI/fabs(Del);
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
      if ((fabs(ans-prev)<epsabs)&&(n>50)) { return ans/M_PI; }
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


    while (Del < 30.) {
      f = lev_int(xbar,Del);
      fout << scientific << Del << "    " << f << endl;
      Del += .005;
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
      f_pos = lev_int(xbar,+Del/mD_over_T);
      f_neg = lev_int(xbar,-Del/mD_over_T);
      f_app = (2./3.)*xbar*exp(-2.*c0*xbar)*pow( fabs(Del), -1.+(4./3.)*xbar/mD_over_T );
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
    double C = 2*(c1+c2) + 1. - GAMMA_E;

    fout << "# columns: Delta/mD/xbar-log(xbar)-C, xbar*mD*f(xbar,Delta)" << endl;

    double f;
    double Del_over_xbar = .01;
    Del_over_xbar = -5.0 + log(xbar)+C;


    while (Del_over_xbar-log(xbar)-C<15) {
      f = lev_int(xbar,Del_over_xbar*xbar)*xbar;
      fout << scientific << Del_over_xbar-log(xbar)-C << "    " << f << endl;
      Del_over_xbar += .01;
    }

    fout.close();
    cout << " finished... " << endl;//*/
                                  //
  }


}

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
double l2f(double x) { return Li2(-exp(-x)) ; }
double l3f(double x) { return Li3(-exp(-x)) ; }
double l1b(double x) { return log(1-exp(-x)); }
double l2b(double x) { return Li2(exp(-x)) ; }
double l3b(double x) { return Li3(exp(-x)) ; }

double R_hard(double eps, double T, double mu, double nf) {
  double mD2 = sqr(mu)*nf/(2.*sqr(M_PI)) + sqr(T)*(1.+nf/6.);
  double Phi_f, Phi_b;

  Phi_b = eps*sqr(T)*( sqr(M_PI)/6. - l2b(eps/T) ) + 2.*cube(T)*( gsl_sf_zeta(3.)-l3b(eps/T) );

  Phi_f = eps*sqr(T)*( l2f((eps+mu)/T)-l2f(+mu/T) ) + 2.*cube(T)*( l3f((eps+mu)/T)-l3f(+mu/T) );
  Phi_f += eps*sqr(T)*( l2f((eps-mu)/T)-l2f(-mu/T) ) + 2.*cube(T)*( l3f((eps-mu)/T)-l3f(-mu/T) );

  double res = 6.*Phi_b + nf*Phi_f;
  res /= cube(eps);
  res /= sqr(2.*M_PI);
  res /= mD2;
  res *= 2*sqr(eps);
  return res;
}


/*--------------------------------------------------------------------*/

int main() {
  using namespace itp;

  init(0.);

  cout << " R(eps=0.1,T=1,nf=3) = " << R_hard(.1,1.,0.,3.) << endl;
  cout << " R(eps=10,T=1,nf=3) = " << R_hard(10.,1.,0.,3.) << endl;
  cout << " R(eps=100,T=1,nf=3) = " << R_hard(100.,1.,0.,3.) << endl;

  /*
  scan_Del(.2,"data/f_T_0_x_0p2.dat");
  scan_Del(.4,"data/f_T_0_x_0p4.dat");
  scan_Del(.6,"data/f_T_0_x_0p6.dat");
  scan_Del(.8,"data/f_T_0_x_0p8.dat");
  scan_Del(1.,"data/f_T_0_x_1p0.dat");
  scan_Del(2.,"data/f_T_0_x_2p0.dat");
  scan_Del(4.,"data/f_T_0_x_4p0.dat");
  scan_Del(8.,"data/f_T_0_x_8p0.dat");//*/


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
