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

  double Delta_min, Delta_max;
  int Delta_N = 100;

  if (T<.01) { Delta_min = .01; } // "effectively" zero temperature case
  else { Delta_min = -5.; } 
  Delta_max = 3e1;
  double Delta_step = (Delta_max-Delta_min)/((double)Delta_N-1);

  //double Delta = Delta_min;

  #pragma omp parallel for
  for (int i=0; i<Delta_N; i++) {
    double Del_tmp = Delta_min + i*Delta_step;
    double f_tmp = lev_int(xbar,Del_tmp);
    if (fabs(Del_tmp)>1e-4) {
      fout << scientific << Del_tmp << "    " << f_tmp << endl;
    }
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

  double Delta_min, Delta_max;
  int Delta_N = 100;

  if (T<.01) { Delta_min = .01; } // "effectively" zero temperature case
  else { Delta_min = -5.0 + log(xbar) + kapp; }
  Delta_max = 3e1;

  double Delta_step = (Delta_max-Delta_min)/((double)Delta_N-1);

  //double f;
  //double Del_over_xbar = .01;
  //if (T>.01) { Del_over_xbar = -5.0 + log(xbar) + kapp; }

  #pragma omp parallel for
  for (int i=0; i<Delta_N; i++) {
    double Del_over_xbar = Delta_min + i*Delta_step;
    double f_tmp = lev_int(xbar,Del_over_xbar*xbar)*xbar;
    if (fabs(Del_over_xbar)>1e-4) {
    	fout << scientific << Del_over_xbar-log(xbar)-kapp << "    " << f_tmp << endl;
    }
  }

  fout.close();
  cout << " finished... " << endl;//*/
}


/*--------------------------------------------------------------------*/

int main() {

  init(.0,1.,3.);

  if (!askToProceed()) return 0;

  cout << endl;

  scan_Del(.1,"data/f_T_0_mu_1_x_0p1.dat");
  scan_Del(.2,"data/f_T_0_mu_1_x_0p2.dat");
  scan_Del(.3,"data/f_T_0_mu_1_x_0p3.dat");
  scan_Del(.4,"data/f_T_0_mu_1_x_0p4.dat");
  scan_Del(.6,"data/f_T_0_mu_1_x_0p6.dat");
  scan_Del(.8,"data/f_T_0_mu_1_x_0p8.dat");
  scan_Del(1, "data/f_T_0_mu_1_x_1p0.dat");
  scan_Del(3, "data/f_T_0_mu_1_x_3p0.dat");
  scan_Del(4, "data/f_T_0_mu_1_x_4p0.dat");
  scan_Del(6, "data/f_T_0_mu_1_x_6p0.dat");
  scan_Del(8, "data/f_T_0_mu_1_x_8p0.dat");
  scan_Del(16,"data/f_T_0_mu_1_x_16p0.dat");

  
  scan_scaled(.1,"data/phi_T_0_mu_1_x_0p1.dat");
  scan_scaled(.2,"data/phi_T_0_mu_1_x_0p2.dat");
  scan_scaled(.3,"data/phi_T_0_mu_1_x_0p3.dat");
  scan_scaled(.4,"data/phi_T_0_mu_1_x_0p4.dat");
  scan_scaled(.6,"data/phi_T_0_mu_1_x_0p6.dat");
  scan_scaled(.8,"data/phi_T_0_mu_1_x_0p8.dat");
  scan_scaled(1, "data/phi_T_0_mu_1_x_1p0.dat");
  scan_scaled(3, "data/phi_T_0_mu_1_x_3p0.dat");
  scan_scaled(4, "data/phi_T_0_mu_1_x_4p0.dat");
  scan_scaled(6, "data/phi_T_0_mu_1_x_6p0.dat");
  scan_scaled(8, "data/phi_T_0_mu_1_x_8p0.dat");
  scan_scaled(16,"data/phi_T_0_mu_1_x_16p0.dat");//*/


  return 0;
}
