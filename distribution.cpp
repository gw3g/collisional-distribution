#include "inc/macros.hpp"
#include "inc/integrators.hpp"
#include "inc/htl.hpp"

using namespace std;

double const c1  =  0.143575;
double const c2  = -0.123714;
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
  static gsl_spline *spline_chi = nullptr;

  void init(double mD_over_T) {
    vd eta_list, hT_list, hL_list, gT_list, gL_list, chi_list;
    string line;

    //ifstream file("data/table_g_h.dat");
    fin.open("data/table_g_h.dat");

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

    fin.close();

    fin.open("data/table_chi.dat");

    while (getline(fin, line)) {
      if (line.empty() || line[0]=='#') continue;
      stringstream ss(line);
      double _eta, _chi1, _chi2, _chi3, _chi4;
      ss >> _eta >> _chi1 >> _chi2 >> _chi3 >> _chi4;
      //eta_list.push_back(_eta); (already done)
      chi_list.push_back(_chi2);
    }

    fin.close();

    auto mk_spline = [](vd &x, vd &y) {
      gsl_spline *spl = gsl_spline_alloc(gsl_interp_cspline, x.sz);
      gsl_spline_init(spl,x.ar,y.ar,x.sz);
      return spl;
    };

    accelerator = gsl_interp_accel_alloc();
    spline_h_T  = mk_spline(eta_list,hT_list);
    spline_h_L  = mk_spline(eta_list,hL_list);
    spline_g_T  = mk_spline(eta_list,gT_list);
    spline_g_L  = mk_spline(eta_list,gL_list);

    spline_chi  = mk_spline(eta_list,chi_list);
  }

  void free() {
    gsl_spline_free(spline_h_T);
    gsl_spline_free(spline_h_L);
    gsl_spline_free(spline_g_T);
    gsl_spline_free(spline_g_L);
    gsl_spline_free(spline_chi);
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
  double chi(double eta) {
    if (eta<1e-3) { return (4.7866666)*sqr(eta); }
    else if (eta>1e3) { 
      double mD_over_T = .3;
      return 4.*( (1./(3.*mD_over_T))*( log(eta) + GAMMA_E ) + (-0.046050746) );
    }
    else { return gsl_spline_eval(spline_chi, eta, accelerator); }
  }


  double _eta_integrand(double eta, void *params) {
    double x = ((double *)params)[0];
    double Del = ((double *)params)[1];
    double h = h_T(eta) + h_L(eta);
    double g = g_T(eta) + g_L(eta);
    double ch = chi(eta);
    //double res = 
    //cos( eta*Del - x*h )*exp( -x*(g-2.*c0) ) - cos( eta*Del );
    //return res*exp( -x*2.*c0 ) ;
    return cos( eta*Del - x*h )*exp( -x*(g+ch) );
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
    Levin series(200,0.);
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
      Del += .05;
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

int main() {
  using namespace itp;

  init(0.3);

  double eta = .9e-3;
  //
  //cout << " eta = " << eta      << " :\n";
  //cout << " hT  = " << h_T(eta) << endl;
  //cout << " hL  = " << h_L(eta) << endl;
  //cout << " gT  = " << g_T(eta) << endl;
  //cout << " gL  = " << g_L(eta) << endl;
  //cout << " chi = " << chi(eta) << endl;

  scan_Del(.2,"data/f_kappa_0p3_x_0p2.dat");
  scan_Del(.4,"data/f_kappa_0p3_x_0p4.dat");
  scan_Del(.6,"data/f_kappa_0p3_x_0p6.dat");
  scan_Del(.8,"data/f_kappa_0p3_x_0p8.dat");
  scan_Del(2.,"data/f_kappa_0p3_x_2p0.dat");
  scan_Del(8.,"data/f_kappa_0p3_x_8p0.dat");

  scan_scaled(.2,"data/scaled_f_kappa_0p3_x_0p2.dat");
  scan_scaled(.4,"data/scaled_f_kappa_0p3_x_0p4.dat");
  scan_scaled(.6,"data/scaled_f_kappa_0p3_x_0p6.dat");
  scan_scaled(.8,"data/scaled_f_kappa_0p3_x_0p8.dat");
  scan_scaled(2.,"data/scaled_f_kappa_0p3_x_2p0.dat");
  scan_scaled(8.,"data/scaled_f_kappa_0p3_x_8p0.dat");


  return 0;
}
