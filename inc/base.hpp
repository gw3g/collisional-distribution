
double S_full(double,double,double,double);
double R_hard(double,double,double,double);

double _DeltaG(double,double,double,double);
double _H(double,double,double,double);
double _G(double,double,double,double);

double _kappa(double,double,double);

double G_large_eta_const(double,double,double);
double DeltaG_small_eta_const(double,double,double);
double H_small_eta_const(double,double,double);

//void tabulate_G_and_H(int,double,double,double);

namespace itp {
  extern double T, mu, nf;
  void init(double,double,double);
  void free();
  double G_func(double);
  double DelG_func(double);
  double H_func(double);
  double lev_int(double,double);
}
