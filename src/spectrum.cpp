#include <iostream>
#include <fstream>
#include <math.h>
#include "base.hpp"

using namespace std;

ofstream fout2;

// file writing (i.t.o. rescaled variables)
void R_fixed_T_mu(double T, double mu, string fname) {

  fout2.open(fname);
  fout2.precision(8);

  fout2 << "# columns: eps/max(pi.T,mu), R(nf=0), R(nf=2), R(nf=3)" << endl;

  double calR_1, calR_2, calR_3;
  double eps = 1e-3;
  double Lambda = max(M_PI*T,mu);

  while (eps<1e3) {
    calR_1 = R_hard(eps*Lambda,T,mu,0);
    calR_2 = R_hard(eps*Lambda,T,mu,2);
    calR_3 = R_hard(eps*Lambda,T,mu,3);
    fout2 << scientific << eps    << "    " 
                       << calR_1 << "    " 
                       << calR_2 << "    " 
                       << calR_3 << endl;
    eps*=1.01;
  }

  fout2.close();
  cout << " finished... " << endl;//*/
                                  //
}

// file writing (i.t.o. rescaled variables)
void scan_S_full(double T, double mu, string fname) {

  fout2.open(fname);
  fout2.precision(8);

  fout2 << "# columns: eps/mD, S[nf=3]" << endl;

  double S_;
  double eps = 1e-2;

  while (eps<1e4) {
    S_ = S_full(eps,T,mu,3);
    fout2 << scientific << eps  << "    " 
                       << S_   << endl;
    eps*=1.01;
  }

  fout2.close();
  cout << " finished... " << endl;//*/
                                  //
}



/*--------------------------------------------------------------------*/

int main() {

  R_fixed_T_mu(3.,3.,"data/calR_T3_mu3.dat"); //*/

}
