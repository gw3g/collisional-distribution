#pragma once
/*--------------------------------------------------------------------*/
// HTL self energies: (re & im parts, z=omega/q)

double Pi_T(double y, double omy) {
  double sty = sqrt(abs(y));
  double ooy = 1./y;
  double res = .0;
  if      (y > .0) res += atanh(sty)/sty;
  else if (y < .0) res += atan(sty)/sty;
  else if (y== .0) res += 1.;
  res = .5*( ooy - (omy/y)*res );
  return res;
}

double dPi_T(double y, double omy) {
  double sty = sqrt(abs(y));
  double ooy = 1./y;
  double res = .0;
  if      (y > .0) res += atanh(sty)/sty;
  else if (y < .0) res += atan(sty)/sty;
  else if (y== .0) res += 1.;
  res = - 3./sqr(2.*y) + (2.+omy)*res/sqr(2.*y);
  return res;
}

double ddPi_T(double y, double omy) {
  double sty = sqrt(abs(y));
  double ooy = 1./y;
  double res = .0;
  if      (y > .0) res += atanh(sty)/sty;
  else if (y < .0) res += atan(sty)/sty;
  else if (y== .0) res += 1.;
  res = - ( (13.*y - 15.) + 3.*(5.-6.*y+sqr(y))*res )/8./sqr(y)/y/omy;
  return res;
}


double Pi_L(double y, double omy) {
  double sty = sqrt(abs(y));
  double ooy = 1./y;
  double res = .0;
  if      (y > .0) res += atanh(sty)/sty;
  else if (y < .0) res += atan(sty)/sty;
  else if (y== .0) res += 1.;
  res = -(omy/y)*( 1.-res );
  return res;
}

double dPi_L(double y, double omy) {
  return - 2.*dPi_T(y,omy);
}

double ddPi_L(double y, double omy) {
  return - 2.*ddPi_T(y,omy);
}


struct Funcd_T {
  double omega2;
  double operator() (const double omy) {
    double y=1.-omy;
    return omega2 - Pi_T(y,omy)/omy;
  }
  double df( const double omy) {
    double y=1.-omy;
    return -(- Pi_T(y,omy) - omy*dPi_T(y,omy) )/sqr(omy);
  }
  double ddf( const double omy) {
    double y=1.-omy;
    return - (2.*Pi_T(y,omy) + 2.*omy*dPi_T(y,omy) + sqr(omy)*ddPi_T(y,omy) )/sqr(omy)/omy;
  }
};

struct Funcd_L {
  double omega2;
  double operator() (const double omy) {
    double y=1.-omy;
    return omega2 - Pi_L(y,omy)/omy;
  }
  double df( const double omy) {
    double y=1.-omy;
    return -(- Pi_L(y,omy) - omy*dPi_L(y,omy) )/sqr(omy);
  }
  double ddf( const double omy) {
    double y=1.-omy;
    return - (2.*Pi_L(y,omy) + 2.*omy*dPi_L(y,omy) + sqr(omy)*ddPi_L(y,omy) )/sqr(omy)/omy;
  }
};

/*--------------------------------------------------------------------*/
// root finder method

template <class T>
double rtsafe(T &funcd, const double x1, const double x2, const double xacc) {
  const int MAXIT=100;
  double xh,xl;
  double fl=funcd(x1);
  double fh=funcd(x2);
  if ((fl > .0 && fh > .0)||(fl < .0 && fh < .0))
    throw("root must be bracketed");
  if (fl == .0) return x1;
  if (fh == .0) return x2;
  if (fl < .0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  double rts=.5*(x1+x2);
  double dxold=abs(x2-x1);
  double dx=dxold;
  double f=funcd(rts);
  double df=funcd.df(rts);
  double ddf=funcd.ddf(rts);
  loop(j,0,MAXIT) {
  //for (int j=0;j<MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > .0) || (abs(2.*f) > abs(dxold*df))) {
      dxold=dx;
      dx=.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts; 
    } else {
      dxold=dx;
      dx=f/df;
      double temp=rts;
      rts -= dx/std::max(.8,std::min(1.2,1.-f*ddf/2./sqr(df)));
      if (temp == rts) return rts;
    }
    if (abs(dx) < xacc) return rts;
    f=funcd(rts);
    df=funcd.df(rts);
    if (f < .0)
      xl=rts;
    else
      xh=rts;
  }
  throw("max interactions exceeded");
}

/*--------------------------------------------------------------------*/
// q2_T(eps) and q2_L(eps)

const double tol=1e-7; // tolerance in root-finder

const double coeffs_T[8] = { // for the small eps expansion of q_T
  1., -4/3./M_PI, -32/9./pow(M_PI,2), 8*(-160 + 9*pow(M_PI,2))/81./pow(M_PI,3),
  4*(-960 + 80*pow(M_PI,2) + pow(M_PI,3))/45./pow(M_PI,4),
 -4*(465920 - 50400*pow(M_PI,2) - 756*pow(M_PI,3) + 729*pow(M_PI,4))/3645./pow(M_PI,5),
 -4*(26808320 - 3548160*pow(M_PI,2) - 57024*pow(M_PI,3) + 95904*pow(M_PI,4) + 729*pow(M_PI,5))/32805./pow(M_PI,6),
 16*(-430080 + 67200*pow(M_PI,2) + 1120*pow(M_PI,3) - 2632*pow(M_PI,4) - 35*pow(M_PI,5) + 15*pow(M_PI,6))/315./pow(M_PI,7)
};

const double coeffs_L[8] = { // for the small eps expansion of q_L
  1., -M_PI/4., 0.5 - 3*pow(M_PI,2)/32., -M_PI*(pow(M_PI,2)-8)/16.,
  ( - 4864 + 3360*pow(M_PI,2) - 315*pow(M_PI,4) )/6144.,
  - 1/15. - 7*M_PI/4. + 5*pow(M_PI,3)/8. - 3*pow(M_PI,5)/64.,
  479/240. - 7*M_PI/60. - 777*pow(M_PI,2)/256. + 3003*pow(M_PI,4)/4096. - 3003*pow(M_PI,6)/65536.,
  (320 + 6592*M_PI - 160*pow(M_PI,2) - 4600*pow(M_PI,3) + 840*pow(M_PI,5) - 45*pow(M_PI,7))/960.
};

double q2_T(double eps) {
  if (eps > 16.) { // use large eps expansion...
    double e2 = sqr(eps), e4 = sqr(e2), e6 = e2*e4, e8 = sqr(e4);
    double L = log(8.*e2);
    double series = .25/e2 
                  + (5.-2.*L)/32./e4 
                  + (27.-2.*(L-9.)*L)/128./e6 
                  + (715.-8.*L*(81.-(L-19.)*L))/2048./e8;
    return sqr(eps*(1.-series));
  } else if (eps < .006) { // use small eps expansion...
    double rho = pow( 4.*sqr(eps)/M_PI, 1/3. );
    double res = 0.;
    loop(i,0,8) {
      res += coeffs_T[i]*pow(rho,i-1.);
    }
    return - sqr(eps*res);
  } else { // otherwise, use root finder:
    Funcd_T fT;
    fT.omega2 = sqr(eps);
    double yT;
    double y1, y2; // -> bracket root
    y1 = -1. + sqrt( 1. + .5/sqr(eps) );
    y2 = .5/sqr(eps);
    //cout << " solving T: " << eps << endl;
    //cout << " y1 = " << y1 << " , y2 = " << y2 << endl;
    yT =  rtsafe(fT,y1,y2,tol);
    double eps2 = eps*eps;
    double qT2 = (1.-yT)*eps2;
    return qT2;
    //double emq2 = pow(eps2-qT2,2.);
    //return 1. - 3./(3.*eps2-1.) - 2.*emq2/(3.*emq2-eps2);
  }
}

double q2_L(double eps) {
  if (eps > 4.) { // use large eps expansion...
    return sqr(eps)*(1.-4.*exp(-2.*(1.+sqr(eps)))); 
  } else if (eps < .02) { // use small eps expansion...
    double res = 0.;
    loop(i,0,8) {
      res += coeffs_L[i]*pow(eps,i);
    }
    return - sqr(res);
  } else { // otherwise, use root finder:
    Funcd_L fL;
    fL.omega2 = sqr(eps);
    double yL;
    double y1, y2; // -> bracket root
    y1 = 4.*exp( -2.*(1.+sqr(eps)) -.1 );
    if (eps<1.) {
      y2 = 1./sqr(eps);
    } else {
      y2 = exp( 2.*(1.-sqr(eps) ) );
    }
    //y2 = 1. - .25/sqr(eps);
    //cout << " solving L: " << eps << endl;
    //cout << " y1 = " << y1 << " , y2 = " << y2 << endl;
    yL =  rtsafe(fL,y1,y2,tol);
    double eps2 = eps*eps;
    double qL2 = (1.-yL)*eps2;
    return qL2;
    //double emq = eps2-qL2;
    //return -1. + 3./(3.*eps2-1.) + 2.*emq/(3.*emq-1.);
  }
}

/*--------------------------------------------------------------------*/

double eval_S_T(double eps, void *params) {
  if (eps<1e-4) { return 1./3. 
                         - (8./9.)*pow(2.*eps/sqr(M_PI),2./3.)
                         - (160./27.)*pow(2.*eps/sqr(M_PI),4./3.)
                         + ( 3.*pow(M_PI,4.)/4. + sqr(4.*M_PI/3.) - 128./3. )*4.*sqr(eps)/pow(M_PI,4.)
                           ; }
  double e2 = sqr(eps), q2 = q2_T(eps);
  return 1. - 3.*e2/(3.*e2-1) - 2.*sqr(e2-q2)/(3.*sqr(e2-q2)-e2);
}

double eval_S_L(double eps, void *params) {
  if (eps<1e-4) { return M_PI*eps/4. 
                         + ( sqr(M_PI/2.) - 4. )*sqr(eps) 
                         + (5.*M_PI*(7.*sqr(M_PI)-48.)/128.)*pow(eps,3.)
    ; }
  double e2 = sqr(eps), q2 = q2_L(eps);
  return - 1. + 3.*e2/(3.*e2-1) + 2.*(e2-q2)/(3.*(e2-q2)-1.);
}

double eval_S_tot(double eps, void *params) {
  //cout << " eps = " << eps << endl;
  if (eps<1e-8) { return 1./3. - (8./9.)*pow(2.*eps/sqr(M_PI),2./3.)+M_PI*eps/4.; }
  double e2 = sqr(eps), q2T = q2_T(eps), q2L = q2_L(eps);
  return  (2/3.)*( 1./(3.*(e2-q2L)-1.) - e2/(3.*sqr(e2-q2T)-e2) );
}

