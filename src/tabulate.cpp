#include <iostream>
#include <fstream>
#include <math.h>
#include "macros.hpp"
#include "base.hpp"
#include <omp.h>

using namespace std;

ofstream fout;

/*--------------------------------------------------------------------*/

#include <gsl/gsl_sf_gamma.h>

void print_asympt(double T, double mu, double nf) {
  double a = DeltaG_small_eta_const(T,mu,nf);
  double b = .5*G_large_eta_const(T,mu,nf);
  double small_c = (16./9.)*pow(2./sqr(M_PI),2./3.)*cos(5*M_PI/6.)*gsl_sf_gamma(5./3.);
  double d = H_small_eta_const(T,mu,nf);
  cout << "\n --> small-eta expansion: \n"  << endl;
  cout << "   G(eta) ~ pi*eta/4"    << endl;
  cout << "  DG(eta) ~ a*eta^2" << endl;
  cout << "        a ≈ " << a << endl;
  cout << "   H(eta) ~ eta/2*[ log(1/eta) + 4*d - Gamma_E + 1 ]" << endl;
  cout << "        d ≈ " << d << endl << endl;
  cout << " --> large-eta expansion: \n"  << endl;
  cout << "   G(eta) ~ 2*b + c*eta^{-5/3}" << endl;
  cout << "        b ≈ " << b << endl;
  cout << "        c = (16/9)*(2/pi^2)^{2/3}*cos(5pi/6)*Gamma(5/3) ≈ " << small_c << endl;
  cout << "  DG(eta) ~ 4/3*T/mD*[ log(eta) + Gamma_E ]" << endl;
  cout << "   H(eta) ~ 2/(3*eta)" << endl;
}

/*--------------------------------------------------------------------*/

void tabulate_G_and_H(int eta_N, double T, double mu, double nf) {

  stringstream filename;
  filename << "data/table_G_H_{T=" << fixed << setprecision(2) << T << ",mu=" << mu << ",nf=" << nf << "}.dat";

  fout.open(filename.str());

  fout.precision(8);
  fout << "# columns:\n";
  fout << "# eta,            G,                DeltaG,           H\n";

  vd eta_list(2*eta_N), G_list(2*eta_N), DelG_list(2*eta_N), H_list(2*eta_N); // lists storage
  //double eta, G_tmp, DeltaG_tmp, H_tmp;
  double eta_min, eta_max;

  //int eta_N = 5000;

  eta_min = 1e-3, eta_max=1e1;
  double delta  = (eta_max-eta_min)/((double)eta_N-1);
  //eta = eta_min;
  cout << "\n --> beginning tabulation (in eta): N = " << eta_N << " , delta_eta = " << delta << endl;
  cout << "                                    T/mD = " << T << " , mu/mD = " << mu << " , nf = " << nf << endl << endl;

  cout << left
       << setw(12) << " eta: "
       << setw(12) << " G: "
       << setw(12) << " Delta G: "
       << setw(12) << " H: " << endl;

  cout << right << fixed << setprecision(5);

  #pragma omp parallel for
  for (int i=0; i<eta_N; i++) {
  //loop(i,0,eta_N) { // first loop is linear scale from eta ~ [0,10]

    //if (eta>1e2) { tolosc = 1e-7; }
    double eta = eta_min + i*delta;  // linear scale
    eta_list[i]  = eta;
    G_list[i] = _G(eta,T,mu,nf);
    DelG_list[i] = _DeltaG(eta,T,mu,nf);
    H_list[i] = _H(eta,T,mu,nf);
    //G_tmp  = _G(eta,T,mu,nf);
    //DeltaG_tmp = _DeltaG(eta,T,mu,nf);
    //H_tmp  = _H(eta,T,mu,nf);

    //eta_list[i]  = eta;
    //G_list[i]    = G_tmp;
    //DelG_list[i] = DeltaG_tmp;
    //H_list[i]    = H_tmp;

//    cout << setw(12) << eta_list[i]
//         << setw(12) << G_list[i]
//         << setw(12) << DelG_list[i]
//         << setw(12) << H_list[i]       << endl;

//    fout << scientific << eta
//         <<     "    " << G_tmp
//         <<     "    " << DeltaG_tmp
//         <<     "    " << H_tmp
//         << endl;

//    eta += delta;
  }

  eta_min = 1e1, eta_max=1e3;
  double ratio  = pow(eta_max/eta_min,1./((double)eta_N-1));
  //eta = eta_min*ratio;

  #pragma omp parallel for
  for (int i=0; i<eta_N-1; i++) {
  //loop(i,1,eta_N) { // second loop is log scale for eta ~ [10,1000]

    //if (eta>1e2) { tolosc = 1e-7; }
    double eta = eta_min * pow(ratio, i+1);  // log scale
    int idx = i + eta_N;                   // store in second half of arrays
    eta_list[idx]  = eta;
    G_list[idx] = _G(eta,T,mu,nf);
    DelG_list[idx] = _DeltaG(eta,T,mu,nf);
    H_list[idx] = _H(eta,T,mu,nf);
    //G_tmp = _G(eta,T,mu,nf);
    //DeltaG_tmp = _DeltaG(eta,T,mu,nf);
    //H_tmp = _H(eta,T,mu,nf);

    //eta_list[i]  = eta;
    //G_list[i]    = G_tmp;
    //DelG_list[i] = DeltaG_tmp;
    //H_list[i]    = H_tmp;

//    cout << setw(12) << eta_list[idx]
//         << setw(12) << G_list[idx]
//         << setw(12) << DelG_list[idx]
//         << setw(12) << H_list[idx]       << endl;

//    fout << scientific << eta
//         <<     "    " << G_tmp
//         <<     "    " << DeltaG_tmp
//         <<     "    " << H_tmp
//         << endl;

//    eta *= ratio;
  }

  loop(i,0,2*eta_N-1) {
    fout << scientific << eta_list[i]
         <<     "    " << G_list[i]
         <<     "    " << DelG_list[i]
         <<     "    " << H_list[i]
         << endl;
  }

  fout.close();
  cout << " finished... " << endl;

}

//*/
/*--------------------------------------------------------------------*/

int main() {

  double T=0.,mu=3,nf=3.;

  print_asympt(T,mu,nf);

  tabulate_G_and_H(5000,T,mu,nf);

  print_asympt(T,5,nf);
  tabulate_G_and_H(5000,T,5,nf);
}
