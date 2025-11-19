#pragma once
/*--------------------------------------------------------------------*/
#include <gsl/gsl_integration.h>
// wrapper for GSL ∫ methods

int    calls = 1e5;

double tolosc=1e-6;
const double epsabs = 1e-7;
const double epsrel = 1e-7;

void integrator(
  double a,
  double b,
  double (*func)(double,void *),
  void (*params),
  double *res, double *err)
{
  /*
   *  generic code: integrate func(x) from x=a to b
   */
  size_t limit = calls;
  gsl_integration_workspace * WS = gsl_integration_workspace_alloc (calls);
  gsl_function f_aux; f_aux.function = func; f_aux.params = params;
  gsl_integration_qag (&f_aux, a, b, epsabs, epsrel, calls, 6, WS, res, err); 
  gsl_integration_workspace_free (WS);
}

void integrate_osc(
  double a,
  double (*func)(double,void *),
  void (*params),
  double *res, double *err,
  enum gsl_integration_qawo_enum sin_or_cos)
{
  /*
   *  generic code: integrate func(x) from x=a to b
   */
  size_t limit = calls;
  gsl_integration_workspace * WS1 = gsl_integration_workspace_alloc (limit);
  gsl_integration_workspace * WS2 = gsl_integration_workspace_alloc (limit);
  double eta = ((double *)params)[0];
  double freq = 2.*M_PI/eta;

  //gsl_integration_qawo_table* qawo_table = gsl_integration_qawo_table_alloc(omega, b - a, GSL_INTEG_SINE, limit);
  gsl_integration_qawo_table* qawo_table = gsl_integration_qawo_table_alloc(eta, freq, sin_or_cos, limit);


  gsl_function f_aux; f_aux.function = func; f_aux.params = params;
  //gsl_integration_qawo(&f_aux, a,    0, 1e-5, limit, WS, qawo_table, res, err);
  gsl_integration_qawf(&f_aux, a,    tolosc, limit, WS1, WS2, qawo_table, res, err);


  gsl_integration_workspace_free (WS1);
  gsl_integration_workspace_free (WS2);
  gsl_integration_qawo_table_free(qawo_table);

}

void integrate_osc_finite(
  double a, 
  double b, 
  double (*func)(double,void *),
  void (*params),
  double *res, double *err,
  enum gsl_integration_qawo_enum sin_or_cos)
{
  /*
   *  generic code: integrate func(x) from x=a to b
   */
  size_t limit = calls;
  gsl_integration_workspace * WS1 = gsl_integration_workspace_alloc (calls);
  double eta = ((double *)params)[0];

  gsl_integration_qawo_table* qawo_table = gsl_integration_qawo_table_alloc(eta, b - a, sin_or_cos, limit);


  gsl_function f_aux; f_aux.function = func; f_aux.params = params;
  gsl_integration_qawo(&f_aux, a,    b, 1e-7, limit, WS1, qawo_table, res, err);
  //gsl_integration_qawf(&f_aux, a,    1e-7, limit, WS1, WS2, qawo_table, res, err);


  gsl_integration_workspace_free (WS1);
  gsl_integration_qawo_table_free(qawo_table);

}

