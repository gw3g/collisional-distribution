#include <iostream>
#include <fstream>
#include <math.h>
#include "macros.hpp"
#include "base.hpp"
#include <omp.h>

using namespace std;
using namespace itp;

string dots = "\033[1;34m:: \033[0m";
string reset = "\033[0m";
string bold  = "\033[1m";

ofstream fout;

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

// file writing (usual parameters)
void scan_Del(double xbar, string fname) {

  fout.open(fname);
  fout.precision(8);
  cout << dots << bold << "Starting scan for {T=" << T << ",mu=" << mu << ",nf=" << nf << "}" << reset << endl;

  fout << "# columns: Delta/mD, mD*f(xbar,Delta)" << endl;

  double f;
  double Del = .01;
  if (T>.01) { Del = -5.; }
  //Del = -1.05013;

  while (Del < 30.) {
  //while (Del < 3.) {
    f = lev_int(xbar,Del);
    fout << scientific << Del << "    " << f << endl;
    Del += .05;
  }

  fout.close();
  cout << " finished... " << endl;//*/
}


// file writing (i.t.o. rescaled variables)
void scan_scaled(double xbar, string fname) {

  fout.open(fname);
  fout.precision(8);
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
}


/*--------------------------------------------------------------------*/

int main() {

  init(.0,10.,3.);

  if (!askToProceed()) return 0;

  cout << endl;

  scan_Del(.2,"data/f_T_0_mu_10_x_0p2.dat");
  scan_Del(.4,"data/f_T_0_mu_10_x_0p4.dat");
  scan_Del(.6,"data/f_T_0_mu_10_x_0p6.dat");
  scan_Del(.8,"data/f_T_0_mu_10_x_0p8.dat");
  scan_Del(4,"data/f_T_0_mu_10_x_4p0.dat");
  scan_Del(16,"data/f_T_0_mu_10_x_16p0.dat");

  return 0;
}
