// [[Rcpp::depends(BH)]]
#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <chrono>
namespace mp = boost::multiprecision;
using mpfr_20 = mp::number<mp::mpfr_float_backend<20,mp::allocate_stack>,mp::et_off>;

/// EXPERIMENT SIMULATION ///

std::vector<int> r_birth_and_death(int n, double beta, double delta, double t){
  int i;
  std::vector<int> res(n);
  NumericVector unif(n);
  unif=runif(n, 0, 1);
  
  double base=1.0+(delta-beta)/(-delta+beta*exp((beta-delta)*t));
  double exponent=(beta*exp((beta-delta)*t)-delta)/((beta-delta)*exp((beta-delta)*t));
  double a=delta*(exp((beta-delta)*t)-1.0)/(beta*exp((beta-delta)*t)-delta);
  
  for (i=0;i<n;++i){
    if (unif[i] > a) {
      res[i]=std::ceil(std::log((1.0-unif[i])*exponent)/log(base));
    } else {
      res[i]=0;
    }
  }
  return res;
}

// derived from rSalvador by Qi Zheng
// [[Rcpp::export]]
Rcpp::List rluria_list(const int n=10, const double rate=1e-8, const double N0=1, const double Nt=1e9, const double mut_fit=1.0,
                       const double death_prob = 0, const double lag=0, const double e=1, const double cv=0, const int trim=0) {
  NumericVector res(n);
  Function gamma("rgamma");
  
  double beta1s=1.0*(1.0-2.0*death_prob)/(1.0-death_prob);
  double beta2, delta2;
  beta2=mut_fit;
  delta2=beta2*death_prob/(1.0-death_prob);
  
  int i,j;
  NumericVector Ttot(n), m(n), mutations(n);
  NumericVector Nts(n);
  if (cv==0) {
    Ttot[0]=log(Nt/N0)/beta1s;
    m[0]=rate*N0*(exp(beta1s*Ttot[0])-1.0)*(1.0-death_prob)/(1.0-2*death_prob);
    Nts[0]=Nt;
    for (i=1;i<n;++i) {
      Ttot[i]=Ttot[0];
      m[i]=m[0];
      Nts[i]=Nt;
    };
  } else {
    Nts=gamma(n, 1/cv/cv, 1/Nt/cv/cv);
    for (i=0;i<n;++i) {
      Ttot[i]=log(Nts[i]/N0)/beta1s;
      m[i]=rate*N0*(exp(beta1s*Ttot[i])-1.0)*(1.0-death_prob)/(1.0-2*death_prob);
    };
  }
  for (i=0;i<n;++i) {
    mutations[i]=rpois(1, m[i])[0];
  }
  
  NumericVector mt(n);
  double epoch, mutants;
  
  for (i=0;i<=n-1;++i) {
    res[i]=0;
  };
  
  for (i=0;i<=n-1;++i) {
    checkUserInterrupt();
    if (mutations[i]==0) {
      res[i]+=0;
    } else {
      mt=runif(mutations[i], 0, 1);
      
      for (j=0; j<mutations[i]; j++) {
        mutants=0;
        
        epoch=log(1.0+mt[j]*(exp(beta1s*Ttot[i])-1.0))/(beta1s);
        double tg=rpois(1, lag)[0]*log(2)/(beta2-delta2);
        
        if (epoch < Ttot[i]-tg) {
          double all=r_birth_and_death(1, beta2, delta2, Ttot[i]-epoch)[0];
          
          mutants+=all;
          
          res[i]+=mutants;
        };
      };
      
      if (e<1) {
        res[i]=rbinom(1, res[i], e)[0];
      }
      if (trim>0) {
        if (res[i]>trim) {
          res[i]=trim;
        };
      };
    };
  };
  
  return Rcpp::List::create(_["mc"] = res, _["nt"] = Nts);
}

// derived from rSalvador by Qi Zheng
// [[Rcpp::export]]
NumericVector rluria_vec(const int n=10, const double rate=1e-8, const double N0=1, const double Nt=1e9, const double mut_fit=1.0,
                         const double death_prob = 0, const double lag=0, const double e=1, const double cv=0, const int trim=0) {
  NumericVector res(n);
  Function gamma("rgamma");
  
  double beta1s=1.0*(1.0-2.0*death_prob)/(1.0-death_prob);
  double beta2, delta2;
  beta2=mut_fit;
  delta2=beta2*death_prob/(1.0-death_prob);
  
  int i,j;
  NumericVector Ttot(n), m(n), mutations(n);
  if (cv==0) {
    Ttot[0]=log(Nt/N0)/beta1s;
    m[0]=rate*N0*(exp(beta1s*Ttot[0])-1.0)*(1.0-death_prob)/(1.0-2*death_prob);
    for (i=1;i<n;++i) {
      Ttot[i]=Ttot[0];
      m[i]=m[0];
    };
  } else {
    NumericVector Nts(n);
    Nts=gamma(n, 1/cv/cv, 1/Nt/cv/cv);
    for (i=0;i<n;++i) {
      Ttot[i]=log(Nts[i]/N0)/beta1s;
      m[i]=rate*N0*(exp(beta1s*Ttot[i])-1.0)*(1.0-death_prob)/(1.0-2*death_prob);
    };
  }
  for (i=0;i<n;++i) {
    mutations[i]=rpois(1, m[i])[0];
  }
  
  NumericVector mt(n);
  double epoch, mutants;
  
  for (i=0;i<=n-1;++i) {
    res[i]=0;
  };
  
  for (i=0;i<=n-1;++i) {
    checkUserInterrupt();
    if (mutations[i]==0) {
      res[i]+=0;
    } else {
      mt=runif(mutations[i], 0, 1);
      
      for (j=0; j<mutations[i]; j++) {
        mutants=0;
        
        epoch=log(1.0+mt[j]*(exp(beta1s*Ttot[i])-1.0))/(beta1s);
        double tg=rpois(1, lag)[0]*log(2)/(beta2-delta2);
        
        if (epoch < Ttot[i]-tg) {
          double all=r_birth_and_death(1, beta2, delta2, Ttot[i]-epoch)[0];
          
          mutants+=all;
          
          res[i]+=mutants;
        };
      };
      
      if (e<1) {
        res[i]=rbinom(1, res[i], e)[0];
      }
      if (trim>0) {
        if (res[i]>trim) {
          res[i]=trim;
        };
      };
    };
  };
  
  return res;
}

/// AULIARY SEQUENCE FOR PROBABILITY COMPUTATION ///

/// GENERAL INTEGRAL FORMULA ///

// [[Rcpp::export]]
std::vector<double> aux_seq_integrate_s(double e=1, double w=1, double d=0, double lag=0, double phi=0, int n=10) {
  using namespace boost::math::quadrature;
  tanh_sinh<double> tanh;
  double r=1./w;
  double l,chi,proxy,factorial;
  int k;
  std::vector<double> res(n+1);
  
  if (e==1 && d==0) {
    auto f = [&r, &k](double x) {
      return (r*pow(x,r)*pow(1.-x,k-1));
    };
    
    res[0]=-exp(-lag+pow(2.,-r)*lag);
    for (k=1;k<=n;++k) {
      proxy=0.0;
      factorial=1.0;
      l=0.0;
      chi=1.0;
      do {
        checkUserInterrupt();
        res[k]=proxy;
        proxy+=factorial*tanh.integrate(f, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
        l++;
        chi*=0.5;
        if ((chi < pow(phi,w)) && (l>1)) {break;}
        factorial*=lag/l;
      } while ((std::abs(res[k]-proxy) > std::numeric_limits<double>::epsilon()) || (l < lag) || (l < 2));
      res[k]*=exp(-lag);
    };
    
  } else {
    auto f1 = [&d, &r, &e](double x) {
      return ((r*pow(x,r-1)*((-1 + e)*x + d*(-e + x)))/(e*(-1 + x) + (-1 + d)*x));
    };
    auto f = [&d, &r, &k, &e](double x) {
      return (r*pow(1.-d,2)*pow(x,r)*e/pow(e+(1.-d-e)*x,2)*pow((1.-x)/(1.+(1.-d-e)/e*x),k-1));
    };
    
    proxy=0.0;
    factorial=1.0;
    l=0.0;
    chi=1.0;
    do {
      checkUserInterrupt();
      res[0]=proxy;
      proxy+=factorial*tanh.integrate(f1, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
      l++;
      chi*=0.5;
      if ((chi < pow(phi,w)) && (l>1)) {break;}
      factorial*=lag/l;
    } while ((std::abs(res[k]-proxy) > std::numeric_limits<double>::epsilon()) || (l < lag) || (l < 2));
    res[0]*=exp(-lag);
    res[0]-=exp(-lag+pow(2.,-r)*lag);
    
    for (k=1;k<=n;++k) {
      proxy=0.0;
      factorial=1.0e0;
      l=0.0;
      chi=1.0;
      do {
        checkUserInterrupt();
        res[k]=proxy;
        proxy+=factorial*tanh.integrate(f, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
        l++;
        chi*=0.5;
        if ((chi < pow(phi,w)) && (l>1)) {break;}
        factorial*=lag/l;
      } while ((std::abs(res[k]-proxy) > std::numeric_limits<double>::epsilon()) || (l < lag) || (l < 2));
      res[k]*=exp(-lag);
    };
  };
  
  return res;
}

/// SEQUENCE FOR LD DISTRIBUTION WITH OPTIONAL PARTIAL PLATING AND DIFFERENTIAL GROWTH ///

// not exported
std::vector<double> aux_seq_exact(double e, double w, int n) {
  std::vector<double> hSeq(n+1);
  int i;
  double r=1.0/w;
  double z = 1.0-e;
  hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    checkUserInterrupt();
    hSeq[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
}

// not exported
std::vector<double> aux_seq_exact_boost(double ee, double w, int n) {
  std::vector<double> hSeq(n+1);
  std::vector<mpfr_20> hSeq2(n+1);
  int i;
  double r=1.0/w;
  double z = 1.0-ee;
  mpfr_20 e=ee;
  hSeq2[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    checkUserInterrupt();
    hSeq2[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  for (i=0; i<=n; ++i) {
    hSeq[i]=hSeq2[i].convert_to<double>();
  };
  return(hSeq);
}

// not exported
std::vector<double> aux_seq_boost(double e, double w, int n) {
  std::vector<double> hSeq(n+1);
  int i;
  double r=1.0/w;
  double z = 1.0-e;
  hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    checkUserInterrupt();
    hSeq[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
}

// not exported
std::vector<mpfr_20> aux_seq_boost(mpfr_20 e, mpfr_20 w, int n) {
  std::vector<mpfr_20> hSeq(n+1);
  int i;
  mpfr_20 r=1.0/w;
  mpfr_20 z = 1.0-e;
  hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(1.0)}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    checkUserInterrupt();
    hSeq[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
}

// partially derived from rSalvador by Qi Zheng
// [[Rcpp::export]]
std::vector<double> aux_seq(double e, double w, int n, int option=0) {
  // NumericVector aux_seq(double e, double w, int n, int option) {
  if (option==1) {return aux_seq_exact(e, w, n);}
  else if (option==2) {return aux_seq_exact_boost(e, w, n);};
  
  if (n<1 || e<=0 || e>1 || w<=0) {return std::vector<double> {0};};
  
  std::vector<double> seq(n+1);
  int i;
  double r=1.0/w;
  double z=1.0-e;
  
  double lowerlimit=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z)*boost::math::beta(1,r+1.0)*pow(e,r)*r;
  double upperlimit=boost::math::hypergeometric_pFq({r, r+1.0}, {r+1.0+n}, z)*boost::math::beta(n,r+1.0)*pow(e,r)*r;
  if (lowerlimit==0 || upperlimit==0) {return std::vector<double> {0};};
  
  if (e==1) {
    seq[0]=-1.0;
    seq[1]=1.0/(1.0+r)*r;
    if (n>1) {
      for (i=2; i<=n; ++i) {
        seq[i]=(i-1.0)/(i+r)*seq[i-1];
      };
    };
  } else if (w!=1 && e!=1) {
    std::vector<double> bSeq(n+1);
    bSeq[0]=-1.0;
    bSeq[1]=1.0/(1.0+r);
    if (n>1) {
      for (i=2; i<=n; ++i) {
        bSeq[i]=(i-1.0)/(i+r)*bSeq[i-1];
      };
    };
    std::vector<double> fSeq(n+1);
    fSeq[0]=0;
    if (n==2) {
      fSeq[1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z);
    } else if (n==3) {
      fSeq[1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z);
      fSeq[2]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+3.0}, z);
    } else if (n>=4) {
      if (e>=0.5) {
        fSeq[n]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+n+1.0}, z);
        fSeq[n-1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+n}, z);
        for (i=(n-1); i>=2; --i) {
          fSeq[i-1]=i*(i+1.0)*z*fSeq[i+1]/(r+i+1.0)/(r+i)/e+(r+i-2.0*i*z)*fSeq[i]/e/(r+i);
        };
      } else {
        fSeq[1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z);
        fSeq[2]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+3.0}, z);
        for (i=2; i<=(n-1); ++i) {
          fSeq[i+1]=(r+1.0+i)/(i*(i+1.0)*z)*(e*(r+i)*fSeq[i-1]-(r+i-2.0*i*z)*fSeq[i]);
        };
      };
    };
    seq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
    for (i=1; i<=n; ++i) {
      seq[i]=pow(e,r)*bSeq[i]*fSeq[i]*r;
    };
  } else if (w==1 && e>=0.5) {
    seq[0]=log(e)*e/z;
    seq[n]=boost::math::hypergeometric_pFq({1.0, 2.0}, {n+2.0}, z)*e/n/(n+1.0);
    if (n>1) {
      for (i=n-1; i>=1; --i) {
        seq[i]=1.0/i/(i+1.0)-seq[i+1]*z/e;
      };
    };
  } else if (w==1 && e<0.5) {
    seq[0]=log(e)*e/z;
    seq[1]=(e/z)*(-1.0-log(e)/z);
    if (n>1) {
      for (i=2; i<=n; ++i) {
        seq[i]=e/z*(1.0/i/(i-1.0)-seq[i-1]);
      };
    };
  } else {
    return(std::vector<double> {0});
  };
  if (abs(seq[1]-lowerlimit)/lowerlimit>1e-9 || abs(seq[n]-upperlimit)/upperlimit>1e-9) {
    return(aux_seq_exact(e, w, n));
  } else {
    return(seq);
  }
}

/// SEQUENCES FOR LD WITH PHENOTYPIC LAG ///

// [[Rcpp::export]]
std::vector<double> aux_seq_lag_s(double e=1, double lag=0, int n=10) {
  int k;
  std::vector<double> hSeq(n+1);
  double max=0;
  
  if (e==1){
    hSeq[0]=-exp(-lag/2);
    for (k=1;k<=n;++k){
      double a=lag;
      double proxy=0;
      double factorial=1;
      int x=0;
      do {
        checkUserInterrupt();
        hSeq[k]=proxy;
        proxy+=factorial*(1.0-pow(1.0-pow(2,-x),k)*(1.0+k*pow(2,-x)))/k/(k+1.0);
        x++;
        if (x>max) {max=x;};
        factorial*=a/x;
      } while (((std::abs(hSeq[k]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
      hSeq[k]*=exp(-a);
    };
  } else {
    int len=n+100;
    std::vector<double> qSeq(len+2);
    double odds=e/(1.0-e);
    double a=lag;
    double L, xi, div;
    // q0
    double proxy=0;
    double factorial=1;
    int x=0;
    do {
      checkUserInterrupt();
      L=pow(2,x);
      xi=e*(1.0-L);
      qSeq[0]=proxy;
      proxy+=factorial*exp(-a)*odds*log(-e*L/(xi-1.0));
      x++;
      if (x>max) {max=x;};
      factorial*=a/x;
    } while (((std::abs(qSeq[0]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
    hSeq[0]=qSeq[0];
    // q1
    if (n >= 1) {
      proxy = 0;
      factorial = 1;
      x = 0;
      do {
        checkUserInterrupt();
        L=pow(2,x);
        xi=e*(1.0-L);
        div=xi/(xi-1.0);
        qSeq[1]=proxy;
        proxy+=factorial*exp(-a)*odds*(1.0/(xi-1.0)-log(-e*L/(xi-1.0)));
        x++;
        if (x>max) {max=x;};
        factorial*=a/x;
      } while (((std::abs(qSeq[1]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
      hSeq[1]=-odds*qSeq[0]+qSeq[1];
    };
    // qn
    if (n >= 2) {
      for (k=2;k<=len;++k){
        proxy = 0;
        factorial = 1;
        x = 0;
        do {
          checkUserInterrupt();
          L=pow(2,x);
          xi=e*(1.0-L);
          div=xi/(xi-1.0);
          qSeq[k]=proxy;
          proxy+=factorial*exp(-a)*odds/k/(k-1)*(1.0-(xi-k)*pow(div,k-1)/(xi-1.0));
          x++;
          if (x>max) {max=x;};
          factorial*=a/x;
        } while (((std::abs(qSeq[k]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
      };
      if (e<0.5) {
        for (k=1; k<=n-1; ++k) {
          hSeq[k+1]=-odds*hSeq[k]+qSeq[k+1];
        };
      } else {
        std::vector<double> seq2(len+1);
        seq2[len]=(1.0-e)*qSeq[len+1];
        for (k=len-1; k>=2; --k) {
          seq2[k]=1/odds*(qSeq[k+1]-seq2[k+1]);
        };
        for (k=2; k<=n; ++k) {
          hSeq[k]=seq2[k];
        };
      };
    };
  };
  
  // Rcout << max << "\n";
  return(hSeq);
  
}

std::vector<mpfr_20> aux_seq_lag_s(mpfr_20 e=1, mpfr_20 lag=0, int n=10) {
  int k;
  std::vector<mpfr_20> hSeq(n+1);
  
  if (e==1){
    hSeq[0]=-exp(-lag/2);
    for (k=1;k<=n;++k){
      mpfr_20 a=lag;
      mpfr_20 proxy=0;
      mpfr_20 factorial=1;
      int x=0;
      do {
        checkUserInterrupt();
        hSeq[k]=proxy;
        proxy+=factorial*exp(-a)*(1.0-pow(1.0-pow(2,-x),k)*(1.0+k*pow(2,-x)))/k/(k+1.0);
        x++;
        factorial*=a/x;
      } while (((mp::abs(hSeq[k]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
    };
  } else {
    int len=n+100;
    std::vector<mpfr_20> qSeq(len+2);
    mpfr_20 odds=e/(1.0-e);
    mpfr_20 a=lag;
    mpfr_20 L, xi, div;
    // q0
    mpfr_20 proxy=0;
    mpfr_20 factorial=1;
    int x=0;
    do {
      checkUserInterrupt();
      L=pow(2,x);
      xi=e*(1.0-L);
      qSeq[0]=proxy;
      proxy+=factorial*exp(-a)*odds*log(-e*L/(xi-1.0));
      x++;
      factorial*=a/x;
    } while (((mp::abs(qSeq[0]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
    hSeq[0]=qSeq[0];
    // q1
    if (n >= 1) {
      proxy = 0;
      factorial = 1;
      x = 0;
      do {
        checkUserInterrupt();
        L=pow(2,x);
        xi=e*(1.0-L);
        div=xi/(xi-1.0);
        qSeq[1]=proxy;
        proxy+=factorial*exp(-a)*odds*(1.0/(xi-1.0)-log(-e*L/(xi-1.0)));
        x++;
        factorial*=a/x;
      } while (((mp::abs(qSeq[1]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
      hSeq[1]=-odds*qSeq[0]+qSeq[1];
    };
    // qn
    if (n >= 2) {
      for (k=2;k<=len;++k){
        proxy = 0;
        factorial = 1;
        x = 0;
        do {
          checkUserInterrupt();
          L=pow(2,x);
          xi=e*(1.0-L);
          div=xi/(xi-1.0);
          qSeq[k]=proxy;
          proxy+=factorial*exp(-a)*odds/k/(k-1)*(1.0-(xi-k)*pow(div,k-1)/(xi-1.0));
          x++;
          factorial*=a/x;
        } while (((mp::abs(qSeq[k]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
      };
      if (e<0.5) {
        for (k=1; k<=n-1; ++k) {
          hSeq[k+1]=-odds*hSeq[k]+qSeq[k+1];
        };
      } else {
        std::vector<mpfr_20> seq2(len+1);
        seq2[len]=(1.0-e)*qSeq[len+1];
        for (k=len-1; k>=2; --k) {
          seq2[k]=1/odds*(qSeq[k+1]-seq2[k+1]);
        };
        for (k=2; k<=n; ++k) {
          hSeq[k]=seq2[k];
        };
      };
    };
  };
  
  return(hSeq);
  
}

// [[Rcpp::export]]
std::vector<double> aux_seq_lag_s_ext(double e, double lag, int n, bool boost=false) {
  std::vector<double> res(n+1);
  if (boost) {
    std::vector<mpfr_20> res2(n+1);
    res2=aux_seq_lag_s(static_cast<mpfr_20>(e), static_cast<mpfr_20>(lag), n);
    for(int k=0;k<=n;++k){
      res[k]=static_cast<double>(res2[k]);
    };
  } else {
    res=aux_seq_lag_s(e, lag, n);
  }
  return res;
}

/// SEQUENCE FOR LD WITH CELL DEATH ///

// not exported
std::vector<double> aux_seq_death(double e, double w, double d, int n) {
  std::vector<double> hSeq(n+1);
  int i;
  double r=1.0/w;
  if (e==1) {
    hSeq[0]=-1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0);
    for (i=1; i<=n; ++i) {
      checkUserInterrupt();
      hSeq[i]=pow(1.0-d,2)*r*boost::math::hypergeometric_pFq({i+1.0, r+1.0}, {r+i+1.0}, d)*boost::math::beta(i,r+1.0);
    };
  } else {
    if (d>=(1.0-2*e)) {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, d)+(1.0-e-d)/e*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, ((e+d-1.0)/e)));
      for (i=1; i<=n; ++i) {
        checkUserInterrupt();
        try {
          hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*boost::math::hypergeometric_pFq({i+1.0, r+1.0}, {i+r+1.0}, ((e+d-1.0)/e));
        } catch (...) {
          Rcpp::stop("Error in aux_seq_death");
        }
      }
    } else {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, d)+(1.0-e-d)/e*(e/(1-d))*boost::math::hypergeometric_pFq({1.0, 1.0}, {r+2.0}, ((e+d-1.0)/(d-1.0))));
      for (i=1; i<=n; ++i) {
        checkUserInterrupt();
        try {
          hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*pow(e/(1.0-d),r+1)*boost::math::hypergeometric_pFq({r+1.0, r}, {i+r+1.0}, ((e+d-1.0)/(d-1.0)));
        } catch (...) {
          Rcpp::stop("Error in aux_seq_death");
        }
      }
    }
  }
  
  return(hSeq);
}

// not exported
std::vector<mpfr_20> aux_seq_death(mpfr_20 e, mpfr_20 w, mpfr_20 d, int n) {
  std::vector<mpfr_20> hSeq(n+1);
  int i;
  mpfr_20 r=1.0/w;
  if (e==1) {
    hSeq[0]=-1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {r+static_cast<mpfr_20>(2.0)}, d)*boost::math::beta(r,2.0);
    for (i=1; i<=n; ++i) {
      checkUserInterrupt();
      hSeq[i]=pow(1.0-d,2)*r*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(i+1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+i+1.0)}, d)*boost::math::beta(i,r+1.0);
    };
  } else {
    if (d>=(1.0-2*e)) {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, d)+(1.0-e-d)/e*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, ((e+d-1.0)/e)));
      for (i=1; i<=n; ++i) {
        checkUserInterrupt();
        hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(i+1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(i+r+1.0)}, ((e+d-1.0)/e));
      }
    } else {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {static_cast<mpfr_20>(r+2.0)}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, d)+(1.0-e-d)/e*(e/(1-d))*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(1.0)}, {static_cast<mpfr_20>(r+2.0)}, ((e+d-1.0)/(d-1.0))));
      for (i=1; i<=n; ++i) {
        checkUserInterrupt();
        hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*pow(e/(1.0-d),r+1)*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(r+1.0), r}, {static_cast<mpfr_20>(i+r+1.0)}, ((e+d-1.0)/(d-1.0)));
      }
    }
  }
  
  return(hSeq);
}

// [[Rcpp::export]]
std::vector<double> aux_seq_death_ext(double e, double w, double d, int n, bool boost=false) {
  std::vector<double> hSeq(n+1);
  int i;
  if (boost==false) {
    std::vector<double> seq(n+1);
    seq=aux_seq_death(e, w, d, n);
    for (i=0;i<=n;++i) {
      hSeq[i]=seq[i];
    };
    return hSeq;
  } else {
    std::vector<mpfr_20> seq(n+1);
    seq=aux_seq_death(static_cast<mpfr_20>(e), static_cast<mpfr_20>(w), static_cast<mpfr_20>(d), n);
    for (i=0;i<=n;++i) {
      hSeq[i]=seq[i].convert_to<double>();
    };
    return hSeq;
  }
}

/// CALCULATION OF PMFS ///

/// PMF FOR LURIA-DELBRUCK APPROXIMATE DISTRIBUTION ///

// derived from rSalvador by Qi Zheng
// [[Rcpp::export]]
std::vector<double> prob_ld(const double m, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<double> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(m*seq[0]); 
  for(n=1;n<=len-1;++n) {
    for(j=1;j<=n;++j){
      prob[n]+=j*seq[j]*prob[n-j];
    };
    prob[n]*=m/n;
  };
  
  return prob;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_ld(const mpfr_20 m, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(m*seq[0]); 
  for(n=1;n<=len-1;++n) {
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob[n]+=j*seq[j]*prob[n-j];
    };
    prob[n]*=m/n;
  };
  
  return prob;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_ld(const mpfr_20 m, std::vector<mpfr_20> &seq, const int len){
  int n,j;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(m*seq[0]); 
  for(n=1;n<=len-1;++n) {
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob[n]+=j*seq[j]*prob[n-j];
    };
    prob[n]*=m/n;
  };
  
  return prob;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<double> prob_ld_deriv(std::vector<double> &seq, std::vector<double> &prob, const int len){
  int n,j;
  std::vector<double> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_ld_deriv(std::vector<mpfr_20> &seq, std::vector<double> &prob, const int len){
  int n,j;
  std::vector<mpfr_20> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_ld_deriv(std::vector<double> &seq, std::vector<mpfr_20> &prob, const int len){
  int n,j;
  std::vector<mpfr_20> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_ld_deriv(std::vector<mpfr_20> &seq, std::vector<mpfr_20> &prob, const int len){
  int n,j;
  std::vector<mpfr_20> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

/// PMF FOR LURIA-DELBRUCK & GAMMA MIXTURE DISTRIBUTION ///

// derived from rSalvador by Qi Zheng
// [[Rcpp::export]]
std::vector<double> xi_seq(const double A, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<double> xi(len);
  
  for(n=1;n<=len-1;++n) {
    xi[n]=0;
  };
  
  xi[0]=log(1.0-A*seq[0]);
  xi[1]=-A*seq[1]/(1.0-A*seq[0]);
  for(n=2;n<=len-1;++n) {
    for(j=1;j<=n-1;++j) {
      xi[n]+=j*xi[j]*seq[n-j];
    };
    xi[n]=-A*(n*seq[n]-xi[n])/(n*(1.0-A*seq[0]));
  };
  
  return xi;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> xi_seq(const mpfr_20 A, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<mpfr_20> xi(len);
  
  for(n=1;n<=len-1;++n) {
    xi[n]=0;
  };
  
  xi[0]=log(1.0-A*seq[0]);
  xi[1]=-A*seq[1]/(1.0-A*seq[0]);
  for(n=2;n<=len-1;++n) {
    checkUserInterrupt();
    for(j=1;j<=n-1;++j) {
      xi[n]+=j*xi[j]*seq[n-j];
    };
    xi[n]=-A*(n*seq[n]-xi[n])/(n*(1.0-A*seq[0]));
  };
  
  return xi;
}

// derived from rSalvador by Qi Zheng
// [[Rcpp::export]]
std::vector<double> prob_b0(const double A, const double k, const double seq0, std::vector<double> &xi, const int len){
  int n,j;
  std::vector<double> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=1.0/pow((1.0-A*seq0),k);
  for(n=1;n<=len-1;++n) {
    for(j=1;j<=n;++j){
      prob[n] += j*xi[j]*prob[n-j];
    };
    prob[n]*=(-k/n);
  };
  
  return prob;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_b0(const mpfr_20 A, const mpfr_20 k, const mpfr_20 seq0, std::vector<mpfr_20> &xi, const int len){
  int n,j;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=1.0/pow((1.0-A*seq0),k);
  for(n=1;n<=len-1;++n) {
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob[n] += j*xi[j]*prob[n-j];
    };
    prob[n]*=(-k/n);
  };
  
  return prob;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<double> prob_b0_deriv1(const double A, const double k, std::vector<double> &seq, std::vector<double> &xi, const int len){
  int n,j;
  std::vector<double> probk1(len),prob1(len);
  double seq0=seq[0];
  
  probk1=prob_b0(A, k+1.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    prob1[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob1[n-1]+=seq[j-1]*probk1[n-j];
    };
  };
  
  return prob1;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_b0_deriv1(const mpfr_20 A, const mpfr_20 k, std::vector<double> &seq, std::vector<mpfr_20> &xi, const int len){
  int n,j;
  std::vector<mpfr_20> probk1(len),prob1(len);
  mpfr_20 seq0=seq[0];
  
  probk1=prob_b0(A, k+1.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    prob1[n]=0;
  };
  
  for(n=1;n<=len;++n){
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob1[n-1]+=seq[j-1]*probk1[n-j];
    };
  };
  
  return prob1;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<double> prob_b0_deriv2(const double A, const double k, std::vector<double> &seq, std::vector<double> &xi, const int len){
  int n,j;
  std::vector<double> probk2(len),conv(len),prob2(len);
  double seq0=seq[0];
  
  probk2=prob_b0(A, k+2.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    conv[n]=0;
    prob2[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      conv[n-1]+=seq[j-1]*seq[n-j];
    };
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob2[n-1]+=conv[j-1]*probk2[n-j];
    };
    prob2[n-1]*=((k+1)/k);
  };
  
  return prob2;
}

// derived from rSalvador by Qi Zheng
// not exported
std::vector<mpfr_20> prob_b0_deriv2(const mpfr_20 A, const mpfr_20 k, std::vector<double> &seq, std::vector<mpfr_20> &xi, const int len){
  int n,j;
  std::vector<mpfr_20> probk2(len),conv(len),prob2(len);
  mpfr_20 seq0=seq[0];
  
  probk2=prob_b0(A, k+2.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    conv[n]=0;
    prob2[n]=0;
  };
  
  for(n=1;n<=len;++n){
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      conv[n-1]+=seq[j-1]*seq[n-j];
    };
  };
  
  for(n=1;n<=len;++n){
    checkUserInterrupt();
    for(j=1;j<=n;++j){
      prob2[n-1]+=conv[j-1]*probk2[n-j];
    };
    prob2[n-1]*=(k+1)/k;
  };
  
  return prob2;
}

/// PROBABILITY FOR POISSON DISTRIBUTION ///
// not exported
std::vector<double> prob_pois(const double m, const int len){
  int n;
  std::vector<double> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(-m); 
  for(n=1;n<=len-1;++n) {
    prob[n]=prob[n-1]*m/n;
  };
  
  return prob;
}

// not exported
std::vector<mpfr_20> prob_pois(const mpfr_20 m, const int len){
  int n;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(-m); 
  for(n=1;n<=len-1;++n) {
    checkUserInterrupt();
    prob[n]=prob[n-1]*m/n;
  };
  
  return prob;
}

/// PMF AND DERIVATIVES FOR FOLD ///

// [[Rcpp::export]]
Rcpp::List calc_probs(const double m, const double cv, std::vector<double> &seq, const int len, const double poisson=0, bool boost=false){
  std::vector<double> prob(len),prob1(len),prob2(len);
  
  if (cv==0) {
    if (boost==false){
      if (poisson==0){
        prob=prob_ld(m, seq, len);
      } else {
        std::vector<double> prob_lc(len), prob_p(len);
        prob_lc=prob_ld(m, seq, len);
        prob_p=prob_pois(poisson, len);
        prob=prob_ld_deriv(prob_lc, prob_p, len);
      };
      prob1=prob_ld_deriv(seq, prob, len);
      prob2=prob_ld_deriv(seq, prob1, len);
    } else {
      std::vector<mpfr_20> xprob(len),xprob1(len),xprob2(len);
      mpfr_20 xm=m;
      int i;
      
      if (poisson==0){
        xprob=prob_ld(xm, seq, len);
      } else {
        std::vector<mpfr_20> xprob_lc(len), xprob_p(len);
        mpfr_20 xpoisson=poisson;
        xprob_lc=prob_ld(xm, seq, len);
        xprob_p=prob_pois(xpoisson, len);
        xprob=prob_ld_deriv(xprob_lc, xprob_p, len);
      };
      xprob1=prob_ld_deriv(seq, xprob, len);
      xprob2=prob_ld_deriv(seq, xprob1, len);
      
      for (i=0;i<=len-1;++i){
        prob[i]=xprob[i].convert_to<double>();
        prob1[i]=xprob1[i].convert_to<double>();
        prob2[i]=xprob2[i].convert_to<double>();
      };
    };
  } else {
    if (boost==false){
      double k=1.0/cv/cv;
      std::vector<double> xi(len);
      double seq0=seq[0];
      double A=m/k;
      xi=xi_seq(A, seq, len);
      if (poisson==0){
        prob=prob_b0(A, k, seq0, xi, len);
        prob1=prob_b0_deriv1(A, k, seq, xi, len);
        prob2=prob_b0_deriv2(A, k, seq, xi, len);
      } else {
        std::vector<double> prob_lc(len), prob_lc_1(len), prob_lc_2(len), prob_p(len);
        prob_lc=prob_b0(A, k, seq0, xi, len);
        prob_p=prob_pois(poisson, len);
        prob=prob_ld_deriv(prob_lc, prob_p, len);
        prob_lc_1=prob_b0_deriv1(A, k, seq, xi, len);
        prob_lc_2=prob_b0_deriv2(A, k, seq, xi, len);
        prob1=prob_ld_deriv(prob_lc_1, prob_p, len);
        prob2=prob_ld_deriv(prob_lc_2, prob_p, len);
      };
    } else {
      std::vector<mpfr_20> xprob(len),xprob1(len),xprob2(len);
      mpfr_20 k=1.0/cv/cv;
      std::vector<mpfr_20> xi(len);
      mpfr_20 seq0=seq[0];
      mpfr_20 A=m/k;
      int i;
      xi=xi_seq(A, seq, len);
      if (poisson==0){
        xprob=prob_b0(A, k, seq0, xi, len);
        xprob1=prob_b0_deriv1(A, k, seq, xi, len);
        xprob2=prob_b0_deriv2(A, k, seq, xi, len);
      } else {
        std::vector<mpfr_20> xprob_lc(len), xprob_lc_1(len), xprob_lc_2(len), xprob_p(len);
        mpfr_20 xpoisson=poisson;
        xprob_lc=prob_b0(A, k, seq0, xi, len);
        xprob_p=prob_pois(xpoisson, len);
        xprob=prob_ld_deriv(xprob_lc, xprob_p, len);
        xprob_lc_1=prob_b0_deriv1(A, k, seq, xi, len);
        xprob_lc_2=prob_b0_deriv2(A, k, seq, xi, len);
        xprob1=prob_ld_deriv(xprob_lc_1, xprob_p, len);
        xprob2=prob_ld_deriv(xprob_lc_2, xprob_p, len);
      };
      for (i=0;i<=len-1;++i){
        prob[i]=xprob[i].convert_to<double>();
        prob1[i]=xprob1[i].convert_to<double>();
        prob2[i]=xprob2[i].convert_to<double>();
      };
    };
  };

  return Rcpp::List::create(Rcpp::Named("prob") = prob,
                            Rcpp::Named("prob1") = prob1,
                            Rcpp::Named("prob2") = prob2);
}

// not exported
std::vector<double> calc_probs_0(const double m, const double cv, std::vector<double> &seq, const int len, const double poisson=0){
  std::vector<double> prob(len);
  
  if (cv==0) {
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<double> prob_lc(len), prob_p(len);
      prob_lc=prob_ld(m, seq, len);
      prob_p=prob_pois(poisson, len);
      prob=prob_ld_deriv(prob_lc, prob_p, len);
    };
  } else {
    double k=1.0/cv/cv;
    std::vector<double> xi(len);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      std::vector<double> prob_lc(len), prob_p(len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob_p=prob_pois(poisson, len);
      prob=prob_ld_deriv(prob_lc, prob_p, len);
    };
  };
  
  return prob;
}

// not exported
std::vector<mpfr_20> calc_probs_0(const mpfr_20 m, const double cv, std::vector<double> &seq, const int len, const double poisson=0){
  std::vector<mpfr_20> prob(len);
  
  if (cv==0) {
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<mpfr_20> prob_lc(len), prob_p(len);
      mpfr_20 xpoisson=poisson;
      prob_lc=prob_ld(m, seq, len);
      prob_p=prob_pois(xpoisson, len);
      prob=prob_ld_deriv(prob_lc, prob_p, len);
    };
  } else {
    mpfr_20 k=1.0/cv/cv;
    std::vector<mpfr_20> xi(len);
    mpfr_20 seq0=seq[0];
    mpfr_20 A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      std::vector<mpfr_20> prob_lc(len), prob_p(len);
      mpfr_20 xpoisson=poisson;
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob_p=prob_pois(xpoisson, len);
      prob=prob_ld_deriv(prob_lc, prob_p, len);
    };
  };
  
  return prob;
}

// [[Rcpp::export]]
double logprob(double m, int len, std::vector<int>& data, std::vector<double>& seq, double k=0, double poisson=0){
  int i,number;
  int samplesize=data.size();
  double sumlogprob=0;
  std::vector<double> prob(len);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<double> xi(len);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    sumlogprob+=log(prob[number]);
  };
  
  return sumlogprob;
}

// [[Rcpp::export]]
double logprob_boost(double xm, int len, std::vector<int>& data, std::vector<double>& seq, double xk=0, double xpoisson=0){
  int i,number;
  int samplesize=data.size();
  mpfr_20 m=xm;
  mpfr_20 sumlogprob=0;
  std::vector<mpfr_20> prob(len);
  
  if (xk==0){
    if (xpoisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      mpfr_20 poisson=xpoisson;
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<mpfr_20> xi(len);
    mpfr_20 seq0=seq[0];
    mpfr_20 k=xk;
    mpfr_20 A=m/k;
    xi=xi_seq(A, seq, len);
    if (xpoisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      mpfr_20 poisson=xpoisson;
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    sumlogprob+=log(prob[number]);
  };
  
  return sumlogprob.convert_to<double>();
}

/// HELPER FUNCTIONS FOR ONE-PARAMETER ESTIMATES ///

// derived from rSalvador by Qi Zheng
void derivatives(double &U, double &J, double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, const bool &fisher=true){
  int i,number;
  int samplesize=data.size();
  std::vector<double> prob(len),prob1(len),prob2(len),subscore(samplesize);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
    prob1=prob_ld_deriv(seq, prob, len);
    if (fisher) {prob2=prob_ld_deriv(seq, prob1, len);};
  } else {
    std::vector<double> xi(len);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
      prob1=prob_b0_deriv1(A, k, seq, xi, len);
      if (fisher) {prob2=prob_b0_deriv2(A, k, seq, xi, len);};
    } else {
      std::vector<double> prob_p(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
      prob_lc_1=prob_b0_deriv1(A, k, seq, xi, len);
      prob1=prob_ld_deriv(prob_p,prob_lc_1,len);
      if (fisher) {
        prob_lc_2=prob_b0_deriv2(A, k, seq, xi, len);
        prob2=prob_ld_deriv(prob_p,prob_lc_2,len);
      };
    };
  };
  
  U=0;
  J=0;
  loglik=0;
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    subscore[i]=(prob1[number]/prob[number]);
    U+=subscore[i];
    if (fisher) {J+=subscore[i]*subscore[i]-prob2[number]/prob[number];};
    loglik+=log(prob[number]);
  };
}

// derived from rSalvador by Qi Zheng
void derivatives_boost(double &U, double &J, double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, const bool &fisher=true){
  int i,number;
  int samplesize=data.size();
  std::vector<mpfr_20> prob(len),prob1(len),prob2(len),subscore(samplesize);
  mpfr_20 xm=m,xk=k,xpoisson=poisson,xU=0,xJ=0,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, len);
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_ld(xm, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
    prob1=prob_ld_deriv(seq, prob, len);
    if (fisher) {prob2=prob_ld_deriv(seq, prob1, len);};
  } else {
    std::vector<mpfr_20> xi(len);
    mpfr_20 seq0=seq[0];
    mpfr_20 xA=xm/xk;
    xi=xi_seq(xA, seq, len);
    if (poisson==0){
      prob=prob_b0(xA, xk, seq0, xi, len);
      prob1=prob_b0_deriv1(xA, xk, seq, xi, len);
      if (fisher) {
        prob2=prob_b0_deriv2(xA, xk, seq, xi, len);
      };
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_b0(xA, xk, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
      prob_lc_1=prob_b0_deriv1(xA, xk, seq, xi, len);
      prob1=prob_ld_deriv(prob_p,prob_lc_1,len);
      if (fisher) {
        prob_lc_2=prob_b0_deriv2(xA, xk, seq, xi, len);
        prob2=prob_ld_deriv(prob_p,prob_lc_2,len);
      };
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    subscore[i]=(prob1[number]/prob[number]);
    xU+=subscore[i];
    if (fisher) {xJ+=subscore[i]*subscore[i]-prob2[number]/prob[number];};
    xloglik+=log(prob[number]);
  };
  
  U=static_cast<double>(xU);
  J=static_cast<double>(xJ);
  loglik=static_cast<double>(xloglik);
}

// derived from rSalvador by Qi Zheng
void loglik(double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson){
  int i,number;
  int samplesize=data.size();
  std::vector<double> prob(len);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<double> xi(len);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  loglik=0;
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    loglik+=log(prob[number]);
  };
}

// derived from rSalvador by Qi Zheng
void loglik_boost(double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson){
  int i,number;
  int samplesize=data.size();
  std::vector<mpfr_20> prob(len);
  mpfr_20 xk=k,xpoisson=poisson,xm=m,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, len);
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_ld(xm, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<mpfr_20> xi(len);
    mpfr_20 seq0=seq[0];
    mpfr_20 xA=xm/xk;
    xi=xi_seq(xA, seq, len);
    if (poisson==0){
      prob=prob_b0(xA, xk, seq0, xi, len);
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_b0(xA, xk, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    xloglik+=log(prob[number]);
  };
  
  loglik=static_cast<double>(xloglik);
}

// inspired by rSalvador by Qi Zheng
// [[Rcpp::export]]
std::vector<double> optim_m(double current_m, double lower_m, double upper_m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, bool verbose=false){
  int iter=0;
  double new_m,new_U,new_J,new_loglik,current_U,current_J,current_loglik,lower_U,lower_J,lower_loglik,upper_U,upper_J,upper_loglik;
  
  derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);
  if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);};
  
  derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);
  if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);};
  
  if(std::abs(upper_U)<1e-6){
    derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);
    if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);};
    return std::vector<double>{upper_m,upper_U,upper_J,upper_loglik};
  };
  
  if(std::abs(lower_U)<1e-6){
    derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);
    if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);};
    return std::vector<double>{lower_m,lower_U,lower_J,lower_loglik};
  };
  
  while(upper_U>0){
    checkUserInterrupt();
    lower_m=upper_m;lower_U=upper_U;lower_J=upper_J;lower_loglik=upper_loglik;
    upper_m*=5;
    derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);
    if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);};
    if(std::abs(upper_U)<1e-6){
      derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);
      if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);};
      return std::vector<double>{upper_m,upper_U,upper_J,upper_loglik};
    };
  };
  while(lower_U<0){
    checkUserInterrupt();
    upper_m=lower_m;upper_U=lower_U;upper_J=lower_J;upper_loglik=lower_loglik;
    lower_m/=10;
    if (lower_m<1e-24) {return {1e-20};};
    derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);
    if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);};
    if(std::abs(lower_U)<1e-6){
      derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);
      if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);};
      return std::vector<double>{lower_m,lower_U,lower_J,lower_loglik};
    };
  };
  
  if(verbose) {Rcout << "boundaries:: m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "boundaries:: U: " << lower_U << " " << upper_U << "\n";};
  if(verbose) {Rcout << "boundaries:: loglik: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m<lower_m) || (current_m>upper_m)) {current_m=(lower_m+upper_m)/2;};
  
  derivatives(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {
    derivatives_boost(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  };
  if(verbose) {Rcout<< "U: " << current_U << " J: " << current_J << " loglik: " << current_loglik << " m: " << current_m << "\n\n";};
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {return {-1.0};};
  
  while(std::abs(current_U)>1.0e-6){
    checkUserInterrupt();
    iter++;
    if (iter>50) {return {-1.0};};
    
    new_m=current_m+current_U/current_J;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    };
    
    derivatives(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {derivatives_boost(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);};
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return {-1.0};};
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if (new_U<0){
      upper_m=new_m;upper_U=new_U;upper_J=new_J;upper_loglik=new_loglik;
    } else if (new_U>0){
      lower_m=new_m;lower_U=new_U;lower_J=new_J;lower_loglik=new_loglik;
    };
    current_m=new_m;current_U=new_U;current_J=new_J;current_loglik=new_loglik;
  };
  
  return std::vector<double>{current_m, current_U, current_J, current_loglik};
}

// inspired by rSalvador by Qi Zheng
// [[Rcpp::export]]
double root_m(double current_m, double lower_m, double upper_m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, const double lalpha, bool verbose=false){
  int iter=0;
  double new_m,new_U,new_J,new_loglik,current_U,current_J,current_loglik,lower_loglik,upper_loglik;
  
  loglik(upper_loglik, upper_m, seq, len, data, k, poisson);
  if(!std::isfinite(upper_loglik)) {loglik_boost(upper_loglik, upper_m, seq, len, data, k, poisson);};
  upper_loglik-=lalpha;
  
  loglik(lower_loglik, lower_m, seq, len, data, k, poisson);
  if(!std::isfinite(lower_loglik)) {loglik_boost(lower_loglik, lower_m, seq, len, data, k, poisson);};
  lower_loglik-=lalpha;
  
  if (std::abs(lower_loglik)<1e-6){return lower_m;};
  if (std::abs(upper_loglik)<1e-6){return upper_m;};
  
  if (upper_loglik*lower_loglik>0){
    if (upper_loglik<lower_loglik){
      int n=0;
      while(upper_loglik>0){
        checkUserInterrupt();
        n++;
        if (n>100) {return -1;}
        lower_m=upper_m;
        lower_loglik=upper_loglik;
        
        upper_m*=5;
        loglik(upper_loglik, upper_m, seq, len, data, k, poisson);
        if(!std::isfinite(upper_loglik)) {loglik_boost(upper_loglik, upper_m, seq, len, data, k, poisson);};
        upper_loglik-=lalpha;
      };
    } else {
      int n=0;
      while(lower_loglik>0){
        checkUserInterrupt();
        n++;
        if (n>100) {return -1;}
        upper_m=lower_m;
        upper_loglik=lower_loglik;
        
        lower_m/=10;
        if (lower_m<1e-24) {return 0;};
        loglik(lower_loglik, lower_m, seq, len, data, k, poisson);
        if(!std::isfinite(lower_loglik)) {loglik_boost(lower_loglik, lower_m, seq, len, data, k, poisson);};
        lower_loglik-=lalpha;
      };
    };
  };
  
  if (std::abs(lower_loglik)<1e-6){return lower_m;};
  if (std::abs(upper_loglik)<1e-6){return upper_m;};
  
  if(verbose) {Rcout << "Boundaries for m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "Boundary log-likelihood: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m>upper_m) || (current_m<lower_m)) {current_m=(lower_m+upper_m)/2;};
  
  derivatives(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {
    derivatives_boost(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  };
  current_loglik-=lalpha;
  if(verbose) {Rcout << "Starting m: " << current_m << " Starting U: " << current_U << " Starting J: " << current_J << " Starting loglik: " << current_loglik << "\n";};
  
  while(std::abs(current_loglik)>1e-6){
    checkUserInterrupt();
    iter++;
    if (iter>50) {return -1.0;};
    
    new_m=current_m-current_loglik/current_U;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    }
    
    derivatives(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {derivatives_boost(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);};
    new_loglik-=lalpha;
    
    if (((lower_loglik<upper_loglik) && (new_loglik<0)) || ((lower_loglik>upper_loglik) && (new_loglik>0))) {
      lower_m=new_m;
      lower_loglik=new_loglik;
    } else if (((lower_loglik<upper_loglik) && (new_loglik>0)) || ((lower_loglik>upper_loglik) && (new_loglik<0))) {
      upper_m=new_m;
      upper_loglik=new_loglik;
    };
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return -1.0;};
    
    current_m=new_m;current_U=new_U;current_loglik=new_loglik;
  };
  return current_m;
}

// inspired by rSalvador by Qi Zheng
// [[Rcpp::export]]
std::vector<double> combo_optim_m(double current_m, double lower_m, double upper_m, const double &R, std::vector<double> &seq1, std::vector<double> &seq2,
                                  const int &len1, const int &len2, std::vector<int> &data1, std::vector<int> &data2,
                                  const double &k1, const double &k2, const double &poisson1, const double &poisson2, bool verbose=false){
  int iter=0;
  double new_m,new_U,new_J,new_loglik,current_U,current_J,current_loglik,lower_U,lower_J,lower_loglik,upper_U,upper_J,upper_loglik,
  U1,J1,loglik1,U2,J2,loglik2;
  
  derivatives(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);};
  derivatives(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);};
  upper_U=U1+R*U2;
  upper_J=J1+R*R*J2;
  upper_loglik=loglik1+loglik2;
  if (std::abs(upper_U)<1e-6) {return std::vector<double>{upper_m,upper_loglik};};
  
  derivatives(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);};
  derivatives(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);};
  lower_U=U1+R*U2;
  lower_J=J1+R*R*J2;
  lower_loglik=loglik1+loglik2;
  if (std::abs(lower_U)<1e-6) {return std::vector<double>{lower_m,lower_loglik};};
  
  while(upper_U>0){
    checkUserInterrupt();
    lower_m=upper_m;lower_U=upper_U;lower_J=upper_J;lower_loglik=upper_loglik;
    upper_m*=5;
    derivatives(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);};
    derivatives(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);};
    upper_U=U1+R*U2;
    upper_J=J1+R*R*J2;
    upper_loglik=loglik1+loglik2;
    if (std::abs(upper_U)<1e-6) {return std::vector<double>{upper_m,upper_loglik};};
  };
  while(lower_U<0){
    checkUserInterrupt();
    upper_m=lower_m;upper_U=lower_U;upper_J=lower_J;upper_loglik=lower_loglik;
    lower_m/=10;
    if (lower_m<1e-24) {return {0};};
    derivatives(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);};
    derivatives(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);};
    lower_U=U1+R*U2;
    lower_J=J1+R*R*J2;
    lower_loglik=loglik1+loglik2;
    if (std::abs(lower_U)<1e-6) {return std::vector<double>{lower_m,lower_loglik};};
  };
  
  if(verbose) {Rcout << "boundaries:: m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "boundaries:: U: " << lower_U << " " << upper_U << "\n";};
  if(verbose) {Rcout << "boundaries:: loglik: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m<lower_m) || (current_m>upper_m)) {current_m=(lower_m+upper_m)/2;};
  
  derivatives(U1, J1, loglik1, current_m, seq1, len1, data1, k1, poisson1, true);
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, current_m, seq1, len1, data1, k1, poisson1, true);};
  derivatives(U2, J2, loglik2, R*current_m, seq2, len2, data2, k2, poisson2, true);
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*current_m, seq2, len2, data2, k2, poisson2, true);};
  current_U=U1+R*U2;
  current_J=J1+R*R*J2;
  current_loglik=loglik1+loglik2;
  
  if(verbose) {Rcout<< "U: " << current_U << " J: " << current_J << " loglik: " << current_loglik << " m: " << current_m << "\n\n";};
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {return {-1.0};};
  
  while(std::abs(current_U)>1.0e-6){
    checkUserInterrupt();
    iter++;
    if (iter>50) {return {-1.0};};
    
    new_m=current_m+current_U/current_J;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    };
    
    derivatives(U1, J1, loglik1, new_m, seq1, len1, data1, k1, poisson1, true);
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, new_m, seq1, len1, data1, k1, poisson1, true);};
    derivatives(U2, J2, loglik2, R*new_m, seq2, len2, data2, k2, poisson2, true);
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*new_m, seq2, len2, data2, k2, poisson2, true);};
    new_U=U1+R*U2;
    new_J=J1+R*R*J2;
    new_loglik=loglik1+loglik2;
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return {-1.0};};
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if (new_U<0){
      upper_m=new_m;upper_U=new_U;upper_J=new_J;upper_loglik=new_loglik;
    } else if (new_U>0){
      lower_m=new_m;lower_U=new_U;lower_J=new_J;lower_loglik=new_loglik;
    };
    current_m=new_m;current_U=new_U;current_J=new_J;current_loglik=new_loglik;
  };
  
  return std::vector<double>{current_m, current_loglik};
}

// [[Rcpp::export]]
std::vector<double> optim_m_every_Nt(double current_m, double lower_m, double upper_m, std::vector<double> &R, Rcpp::List &seqs,
                                     std::vector<int> &data, const double &k, std::vector<int> &poisson, bool verbose=false){
  int iter=0;
  double new_m=0,new_U=0,new_J=0,new_loglik=0,current_U=0,current_J=0,current_loglik=0,
    lower_U=0,lower_J=0,lower_loglik=0,upper_U=0,upper_J=0,upper_loglik=0,
    proxy_U=0,proxy_J=0,proxy_loglik=0,pois=0;
  std::vector <int> subdata(1);
  int size = data.size();
  
  proxy_U=0; proxy_J=0; proxy_loglik=0; upper_U=0; upper_J=0; upper_loglik=0;
  for (int i=0; i<size; ++i) {
    subdata[0] = data[i];
    pois = poisson[i];
    std::vector <double> seq = seqs[i];
    derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois, true);
    if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois, true);};
    upper_U += R[i]*proxy_U;
    upper_J += R[i]*R[i]*proxy_J;
    upper_loglik += proxy_loglik;
  }
  if (std::abs(upper_U)<1e-6) {return std::vector<double>{upper_m,upper_loglik};};
  
  proxy_U=0; proxy_J=0; proxy_loglik=0; lower_U=0; lower_J=0; lower_loglik=0;
  for (int i=0; i<size; ++i) {
    subdata[0] = data[i];
    pois = poisson[i];
    std::vector <double> seq = seqs[i];
    derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois, true);
    if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois, true);};
    lower_U += R[i]*proxy_U;
    lower_J += R[i]*R[i]*proxy_J;
    lower_loglik += proxy_loglik;
  }
  if (std::abs(lower_U)<1e-6) {return std::vector<double>{lower_m,lower_loglik};};
  
  while(upper_U>0){
    checkUserInterrupt();
    lower_m=upper_m;lower_U=upper_U;lower_J=upper_J;lower_loglik=upper_loglik;
    upper_m*=5;
    proxy_U=0; proxy_J=0; proxy_loglik=0; upper_U=0; upper_J=0; upper_loglik=0;
    for (int i=0; i<size; ++i) {
      subdata[0] = data[i];
      pois = poisson[i];
      std::vector <double> seq = seqs[i];
      derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois, true);
      if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois, true);};
      upper_U += R[i]*proxy_U;
      upper_J += R[i]*R[i]*proxy_J;
      upper_loglik += proxy_loglik;
    }
    if (std::abs(upper_U)<1e-6) {return std::vector<double>{upper_m,upper_loglik};};
  };
  while(lower_U<0){
    checkUserInterrupt();
    upper_m=lower_m;upper_U=lower_U;upper_J=lower_J;upper_loglik=lower_loglik;
    lower_m/=10;
    if (lower_m<1e-24) {return {0};};
    proxy_U=0; proxy_J=0; proxy_loglik=0; lower_U=0; lower_J=0; lower_loglik=0;
    for (int i=0; i<size; ++i) {
      subdata[0] = data[i];
      pois = poisson[i];
      std::vector <double> seq = seqs[i];
      derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois, true);
      if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois, true);};
      lower_U += R[i]*proxy_U;
      lower_J += R[i]*R[i]*proxy_J;
      lower_loglik += proxy_loglik;
    }
    if (std::abs(lower_U)<1e-6) {return std::vector<double>{lower_m,lower_loglik};};
  };
  
  if(verbose) {Rcout << "boundaries:: m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "boundaries:: U: " << lower_U << " " << upper_U << "\n";};
  if(verbose) {Rcout << "boundaries:: loglik: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m<lower_m) || (current_m>upper_m)) {current_m=(lower_m+upper_m)/2;};
  
  proxy_U=0; proxy_J=0; proxy_loglik=0;
  for (int i=0; i<size; ++i) {
    subdata[0] = data[i];
    pois = poisson[i];
    std::vector <double> seq = seqs[i];
    derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*current_m, seq, data[i]+1, subdata, k, pois, true);
    if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*current_m, seq, data[i]+1, subdata, k, pois, true);};
    current_U += R[i]*proxy_U;
    current_J += R[i]*R[i]*proxy_J;
    current_loglik += proxy_loglik;
  }
  
  if(verbose) {Rcout<< "U: " << current_U << " J: " << current_J << " loglik: " << current_loglik << " m: " << current_m << "\n\n";};
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {return {-1.0};};
  
  while(std::abs(current_U)>1.0e-6){
    checkUserInterrupt();
    iter++;
    if (iter>70) {return {-1.0};};
    
    new_m=current_m+current_U/current_J;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    };
    
    proxy_U=0; proxy_J=0; proxy_loglik=0; new_U=0; new_J=0; new_loglik=0;
    for (int i=0; i<size; ++i) {
      subdata[0] = data[i];
      pois = poisson[i];
      std::vector <double> seq = seqs[i];
      derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*new_m, seq, data[i]+1, subdata, k, pois, true);
      if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*new_m, seq, data[i]+1, subdata, k, pois, true);};
      new_U += R[i]*proxy_U;
      new_J += R[i]*R[i]*proxy_J;
      new_loglik += proxy_loglik;
    }
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return {-1.0};};
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if (new_U<0){
      upper_m=new_m;upper_U=new_U;upper_J=new_J;upper_loglik=new_loglik;
    } else if (new_U>0){
      lower_m=new_m;lower_U=new_U;lower_J=new_J;lower_loglik=new_loglik;
    };
    current_m=new_m;current_U=new_U;current_J=new_J;current_loglik=new_loglik;
  };
  
  return std::vector<double>{current_m, current_U, current_J, current_loglik};
}

// [[Rcpp::export]]
double root_m_every_Nt(double current_m, double lower_m, double upper_m, Rcpp::List &seqs, std::vector<double> &R,
                       std::vector<int> &data, const double &k, std::vector<int> &poisson, const double lalpha, bool verbose=false){
  int iter=0;
  double new_m=0,new_U=0,new_J=0,new_loglik=0,current_U=0,current_J=0,current_loglik=0,
    lower_loglik=0,upper_loglik=0,proxy_U=0,proxy_J=0,proxy_loglik=0,pois=0;
  std::vector <int> subdata(1);
  int size = data.size();
  
  proxy_loglik=0; upper_loglik=0;
  for (int i=0; i<size; ++i) {
    subdata[0] = data[i];
    pois = poisson[i];
    std::vector <double> seq = seqs[i];
    loglik(proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois);
    if(!std::isfinite(upper_loglik)) {loglik_boost(proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois);};
    upper_loglik += proxy_loglik;
  }
  upper_loglik-=lalpha;
  
  proxy_loglik=0; lower_loglik=0;
  for (int i=0; i<size; ++i) {
    subdata[0] = data[i];
    pois = poisson[i];
    std::vector <double> seq = seqs[i];
    loglik(proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois);
    if(!std::isfinite(lower_loglik)) {loglik_boost(proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois);};
    lower_loglik += proxy_loglik;
  }
  lower_loglik-=lalpha;
  
  if (std::abs(lower_loglik)<1e-6){return lower_m;};
  if (std::abs(upper_loglik)<1e-6){return upper_m;};
  
  if (upper_loglik*lower_loglik>0){
    if (upper_loglik<lower_loglik){
      int n=0;
      while(upper_loglik>0){
        checkUserInterrupt();
        n++;
        if (n>100) {return -1;}
        lower_m=upper_m;
        lower_loglik=upper_loglik;
        
        upper_m*=5;
        proxy_loglik=0; upper_loglik=0;
        for (int i=0; i<size; ++i) {
          subdata[0] = data[i];
          pois = poisson[i];
          std::vector <double> seq = seqs[i];
          loglik(proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois);
          if(!std::isfinite(upper_loglik)) {loglik_boost(proxy_loglik, R[i]*upper_m, seq, data[i]+1, subdata, k, pois);};
          upper_loglik += proxy_loglik;
        }
        upper_loglik-=lalpha;
      };
    } else {
      int n=0;
      while(lower_loglik>0){
        checkUserInterrupt();
        n++;
        if (n>100) {return -1;}
        upper_m=lower_m;
        upper_loglik=lower_loglik;
        
        lower_m/=10;
        if (lower_m<1e-24) {return 0;};
        proxy_loglik=0; lower_loglik=0;
        for (int i=0; i<size; ++i) {
          subdata[0] = data[i];
          pois = poisson[i];
          std::vector <double> seq = seqs[i];
          loglik(proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois);
          if(!std::isfinite(lower_loglik)) {loglik_boost(proxy_loglik, R[i]*lower_m, seq, data[i]+1, subdata, k, pois);};
          lower_loglik += proxy_loglik;
        }
        lower_loglik-=lalpha;
      };
    };
  };
  
  if (std::abs(lower_loglik)<1e-6){return lower_m;};
  if (std::abs(upper_loglik)<1e-6){return upper_m;};
  
  if(verbose) {Rcout << "Boundaries for m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "Boundary log-likelihood: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m>upper_m) || (current_m<lower_m)) {current_m=(lower_m+upper_m)/2;};
  
  proxy_U=0; proxy_J=0; proxy_loglik=0;
  for (int i=0; i<size; ++i) {
    subdata[0] = data[i];
    pois = poisson[i];
    std::vector <double> seq = seqs[i];
    derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*current_m, seq, data[i]+1, subdata, k, pois, true);
    if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*current_m, seq, data[i]+1, subdata, k, pois, true);};
    current_U += R[i]*proxy_U;
    current_J += R[i]*R[i]*proxy_J;
    current_loglik += proxy_loglik;
  }
  current_loglik-=lalpha;
  if(verbose) {Rcout << "Starting m: " << current_m << " Starting U: " << current_U << " Starting J: " << current_J << " Starting loglik: " << current_loglik << "\n";};
  
  while(std::abs(current_loglik)>1e-6){
    checkUserInterrupt();
    iter++;
    if (iter>50) {return -1.0;};
    
    new_m=current_m-current_loglik/current_U;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    }
    
    proxy_U=0; proxy_J=0; proxy_loglik=0; new_U=0; new_J=0; new_loglik=0;
    for (int i=0; i<size; ++i) {
      subdata[0] = data[i];
      pois = poisson[i];
      std::vector <double> seq = seqs[i];
      derivatives(proxy_U, proxy_J, proxy_loglik, R[i]*new_m, seq, data[i]+1, subdata, k, pois, true);
      if(!std::isfinite(proxy_U) || !std::isfinite(proxy_J) || !std::isfinite(proxy_loglik)) {derivatives_boost(proxy_U, proxy_J, proxy_loglik, R[i]*new_m, seq, data[i]+1, subdata, k, pois, true);};
      new_U += R[i]*proxy_U;
      new_J += R[i]*R[i]*proxy_J;
      new_loglik += proxy_loglik;
    }
    new_loglik-=lalpha;
    
    if (((lower_loglik<upper_loglik) && (new_loglik<0)) || ((lower_loglik>upper_loglik) && (new_loglik>0))) {
      lower_m=new_m;
      lower_loglik=new_loglik;
    } else if (((lower_loglik<upper_loglik) && (new_loglik>0)) || ((lower_loglik>upper_loglik) && (new_loglik<0))) {
      upper_m=new_m;
      upper_loglik=new_loglik;
    };
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return -1.0;};
    
    current_m=new_m;current_U=new_U;current_loglik=new_loglik;
  };
  return current_m;
}

/// FUNCTIONS FOR POWER AND SAMPLE SIZE ///

/// GENERAL INTEGRAL FORMULA ///

// [[Rcpp::export]]
std::vector<double> aux_seq_integrate_s_n0(double e=1, double w=1, double d=0, double lag=0, double phi=0, int n=10, int n0=0) {
  // std::vector<double> aux_seq_integrate(double e, double w, double d, double lag, double phi, int n) {
  using namespace boost::math::quadrature;
  tanh_sinh<double> tanh;
  double r=1./w;
  double l,chi,proxy,factorial;
  int k;
  std::vector<double> res(n-n0+1);
  
  if (e==1 && d==0) {
    auto f = [&r, &k](double x) {
      return (r*pow(x,r)*pow(1.-x,k-1));
    };
    
    if (n0==0) {res[0]=-exp(-lag+pow(2.,-r)*lag);};
    for (k=std::max(n0,1);k<=n;++k) {
      proxy=0.0;
      factorial=1.0;
      l=0.0;
      chi=1.0;
      do {
        checkUserInterrupt();
        res[k-n0]=proxy;
        proxy+=factorial*tanh.integrate(f, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
        l++;
        chi*=0.5;
        if ((chi < pow(phi,w)) && (l>1)) {break;}
        factorial*=lag/l;
      } while ((std::abs(res[k]-proxy) > std::numeric_limits<double>::epsilon()) || (l < lag) || (l < 2));
      res[k-n0]*=exp(-lag);
    };
    
  } else {
    auto f1 = [&d, &r, &e](double x) {
      return ((r*pow(x,r-1)*((-1 + e)*x + d*(-e + x)))/(e*(-1 + x) + (-1 + d)*x));
    };
    auto f = [&d, &r, &k, &e](double x) {
      return (r*pow(1.-d,2)*pow(x,r)*e/pow(e+(1.-d-e)*x,2)*pow((1.-x)/(1.+(1.-d-e)/e*x),k-1));
    };
    
    if (n0==0) {
      proxy=0.0;
      factorial=1.0;
      l=0.0;
      chi=1.0;
      do {
        checkUserInterrupt();
        res[0]=proxy;
        proxy+=factorial*tanh.integrate(f1, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
        l++;
        chi*=0.5;
        if ((chi < pow(phi,w)) && (l>1)) {break;}
        factorial*=lag/l;
      } while ((std::abs(res[k]-proxy) > std::numeric_limits<double>::epsilon()) || (l < lag) || (l < 2));
      res[0]*=exp(-lag);
      res[0]-=exp(-lag+pow(2.,-r)*lag);
    };
    
    for (k=std::max(n0,1);k<=n;++k) {
      proxy=0.0;
      factorial=1.0e0;
      l=0.0;
      chi=1.0;
      do {
        checkUserInterrupt();
        res[k-n0]=proxy;
        proxy+=factorial*tanh.integrate(f, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
        l++;
        chi*=0.5;
        if ((chi < pow(phi,w)) && (l>1)) {break;}
        factorial*=lag/l;
      } while ((std::abs(res[k]-proxy) > std::numeric_limits<double>::epsilon()) || (l < lag) || (l < 2));
      res[k-n0]*=exp(-lag);
    };
  };
  return res;
}

/// SEQUENCE FOR LD DISTRIBUTION WITH OPTIONAL PARTIAL PLATING AND DIFFERENTIAL GROWTH ///

// not exported
std::vector<double> aux_seq_n0(double e=1, double w=1, int n=10, int n0=0) {
  std::vector<double> hSeq(n-n0+1);
  int i;
  double r=1.0/w;
  double z = 1.0-e;
  if (n0==0) {hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);};
  for (i=std::max(1,n0); i<=n; ++i) {
    hSeq[i-n0]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
}

// not exported
std::vector<mpfr_20> aux_seq_n0(mpfr_20 e=1, mpfr_20 w=1, int n=10, int n0=0) {
  std::vector<mpfr_20> hSeq(n-n0+1);
  int i;
  mpfr_20 r=1.0/w;
  mpfr_20 z = 1.0-e;
  if (n0==0) {hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(1.0)}, {2.0 + r}, z)/(1.0+r);};
  for (i=std::max(1,n0); i<=n; ++i) {
    hSeq[i-n0]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
  
}

// [[Rcpp::export]]
std::vector<double> aux_seq_ext_n0(double e=1, double w=0, int n=10, bool boost=false, int n0=0) {
  std::vector<double> res(n-n0+1);
  if (boost) {
    std::vector<mpfr_20> res2(n-n0+1);
    res2=aux_seq_n0(static_cast<mpfr_20>(e), static_cast<mpfr_20>(w), n, n0);
    for(int k=0;k<=(n-n0);++k){
      res[k]=static_cast<double>(res2[k]);
    };
  } else {
    res=aux_seq_n0(e, w, n, n0);
  }
  return res;
}

/// SEQUENCES FOR LD WITH PHENOTYPIC LAG ///

// not exported
std::vector<double> aux_seq_lag_s_n0(double e=1, double lag=0, int n=10, int n0=0) {
  int k;
  std::vector<double> hSeq(n-n0+1);
  double max=0;
  
  if (e==1){
    if (n0==0) {hSeq[0]=-exp(-lag/2);};
    for (k=std::max(1,n0);k<=n;++k){
      double a=lag;
      double proxy=0;
      double factorial=1;
      int x=0;
      do {
        hSeq[k-n0]=proxy;
        proxy+=factorial*(1.0-pow(1.0-pow(2,-x),k)*(1.0+k*pow(2,-x)))/k/(k+1.0);
        x++;
        if (x>max) {max=x;};
        factorial*=a/x;
      } while (((std::abs(hSeq[k-n0]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
      hSeq[k-n0]*=exp(-a);
    };
  } else {
    int len=n+100;
    std::vector<double> qSeq(len+2);
    double odds=e/(1.0-e);
    double a=lag;
    double L, xi, div;
    if (n0<=1){
      // q0
      double proxy=0;
      double factorial=1;
      int x=0;
      do {
        L=pow(2,x);
        xi=e*(1.0-L);
        qSeq[0]=proxy;
        proxy+=factorial*exp(-a)*odds*log(-e*L/(xi-1.0));
        x++;
        if (x>max) {max=x;};
        factorial*=a/x;
      } while (((std::abs(qSeq[0]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
    };
    if (n0==0){
      hSeq[0]=qSeq[0];
    };
    // q1
    if (n >= 1) {
      double proxy = 0;
      double factorial = 1;
      int x = 0;
      do {
        L=pow(2,x);
        xi=e*(1.0-L);
        div=xi/(xi-1.0);
        qSeq[1]=proxy;
        proxy+=factorial*exp(-a)*odds*(1.0/(xi-1.0)-log(-e*L/(xi-1.0)));
        x++;
        if (x>max) {max=x;};
        factorial*=a/x;
      } while (((std::abs(qSeq[1]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
    };
    if (n0<=1) {
      hSeq[1-n0]=-odds*qSeq[0]+qSeq[1];
    };
    // qn
    if (n >= 2) {
      for (k=std::max(n0,2);k<=len;++k){
        double proxy = 0;
        double factorial = 1;
        int x = 0;
        do {
          L=pow(2,x);
          xi=e*(1.0-L);
          div=xi/(xi-1.0);
          qSeq[k]=proxy;
          proxy+=factorial*exp(-a)*odds/k/(k-1)*(1.0-(xi-k)*pow(div,k-1)/(xi-1.0));
          x++;
          if (x>max) {max=x;};
          factorial*=a/x;
        } while (((std::abs(qSeq[k]-proxy) > std::numeric_limits<double>::epsilon()) || (x < a)) || (x < 2));
      };
      if (e<0.5) {
        for (k=std::max(n0-1,1); k<=n-1; ++k) {
          hSeq[k+1-n0]=-odds*hSeq[k-n0]+qSeq[k+1];
        };
      } else {
        std::vector<double> seq2(len+1);
        seq2[len]=(1.0-e)*qSeq[len+1];
        for (k=len-1; k>=std::max(n0,2); --k) {
          seq2[k]=1/odds*(qSeq[k+1]-seq2[k+1]);
        };
        for (k=std::max(n0,2); k<=n; ++k) {
          hSeq[k-n0]=seq2[k];
        };
      };
    };
  };
  return(hSeq);
}

// not exported
std::vector<mpfr_20> aux_seq_lag_s_n0(mpfr_20 e=1, mpfr_20 lag=0, int n=10, int n0=0) {
  int k;
  std::vector<mpfr_20> hSeq(n-n0+1);
  
  if (e==1){
    if (n0==0) {hSeq[0]=-exp(-lag/2);};
    for (k=std::max(1,n0);k<=n;++k){
      mpfr_20 a=lag;
      mpfr_20 proxy=0;
      mpfr_20 factorial=1;
      int x=0;
      do {
        hSeq[k-n0]=proxy;
        proxy+=factorial*exp(-a)*(1.0-pow(1.0-pow(2,-x),k)*(1.0+k*pow(2,-x)))/k/(k+1.0);
        x++;
        factorial*=a/x;
      } while (((mp::abs(hSeq[k-n0]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
    };
  } else {
    int len=n+100;
    std::vector<mpfr_20> qSeq(len+2);
    mpfr_20 odds=e/(1.0-e);
    mpfr_20 a=lag;
    mpfr_20 L, xi, div;
    // q0
    mpfr_20 proxy=0;
    mpfr_20 factorial=1;
    int x=0;
    do {
      L=pow(2,x);
      xi=e*(1.0-L);
      qSeq[0]=proxy;
      proxy+=factorial*exp(-a)*odds*log(-e*L/(xi-1.0));
      x++;
      factorial*=a/x;
    } while (((mp::abs(qSeq[0]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
    if (n0==0){
      hSeq[0]=qSeq[0];
    };
    // q1
    if (n >= 1) {
      mpfr_20 proxy = 0;
      mpfr_20 factorial = 1;
      mpfr_20 x = 0;
      do {
        L=pow(2,x);
        xi=e*(1.0-L);
        div=xi/(xi-1.0);
        qSeq[1]=proxy;
        proxy+=factorial*exp(-a)*odds*(1.0/(xi-1.0)-log(-e*L/(xi-1.0)));
        x++;
        factorial*=a/x;
      } while (((mp::abs(qSeq[1]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
    };
    if (n0<=1) {
      hSeq[1-n0]=-odds*qSeq[0]+qSeq[1];
    };
    // qn
    if (n >= 2) {
      for (k=std::max(n0,2);k<=len;++k){
        mpfr_20 proxy = 0;
        mpfr_20 factorial = 1;
        mpfr_20 x = 0;
        do {
          L=pow(2,x);
          xi=e*(1.0-L);
          div=xi/(xi-1.0);
          qSeq[k]=proxy;
          proxy+=factorial*exp(-a)*odds/k/(k-1)*(1.0-(xi-k)*pow(div,k-1)/(xi-1.0));
          x++;
          factorial*=a/x;
        } while (((mp::abs(qSeq[k]-proxy) > std::numeric_limits<mpfr_20>::epsilon()) || (x < a)) || (x < 2));
      };
      if (e<0.5) {
        for (k=std::max(n0-1,1); k<=n-1; ++k) {
          hSeq[k+1-n0]=-odds*hSeq[k-n0]+qSeq[k+1];
        };
      } else {
        std::vector<mpfr_20> seq2(len+1);
        seq2[len]=(1.0-e)*qSeq[len+1];
        for (k=len-1; k>=std::max(n0,2); --k) {
          seq2[k]=1/odds*(qSeq[k+1]-seq2[k+1]);
        };
        for (k=std::max(n0,2); k<=n; ++k) {
          hSeq[k-n0]=seq2[k];
        };
      };
    };
  };
  return(hSeq);
}

// [[Rcpp::export]]
std::vector<double> aux_seq_lag_s_ext_n0(double e=1, double lag=0, int n=10, bool boost=false, int n0=0) {
  std::vector<double> res(n-n0+1);
  if (boost) {
    std::vector<mpfr_20> res2(n-n0+1);
    res2=aux_seq_lag_s_n0(static_cast<mpfr_20>(e), static_cast<mpfr_20>(lag), n, n0);
    for(int k=0;k<=(n-n0);++k){
      res[k]=static_cast<double>(res2[k]);
    };
  } else {
    res=aux_seq_lag_s_n0(e, lag, n, n0);
  }
  return res;
}

/// SEQUENCE FOR LD WITH CELL DEATH ///

// not exported
std::vector<double> aux_seq_death_n0(double e, double w, double d, int n, int n0=0) {
  std::vector<double> hSeq(n+1-n0);
  int i;
  double r=1.0/w;
  if (e==1) {
    if (n0==0) {hSeq[0]=-1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0);};
    for (i=std::max(n0,1); i<=n; ++i) {
      hSeq[i-n0]=pow(1.0-d,2)*r*boost::math::hypergeometric_pFq({i+1.0, r+1.0}, {r+i+1.0}, d)*boost::math::beta(i,r+1.0);
    };
  } else {
    if (d>=(1.0-2*e)) {
      if (n0==0) {hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, d)+(1.0-e-d)/e*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, ((e+d-1.0)/e)));};
      for (i=std::max(n0,1); i<=n; ++i) {
        try {
          hSeq[i-n0] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*boost::math::hypergeometric_pFq({i+1.0, r+1.0}, {i+r+1.0}, ((e+d-1.0)/e));
        } catch (...) {
          Rcpp::stop("Error");
        }
      }
    } else {
      if (n0==0) {hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, d)+(1.0-e-d)/e*(e/(1-d))*boost::math::hypergeometric_pFq({1.0, 1.0}, {r+2.0}, ((e+d-1.0)/(d-1.0))));};
      for (i=std::max(n0,1); i<=n; ++i) {
        try {
          hSeq[i-n0] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*pow(e/(1.0-d),r+1)*boost::math::hypergeometric_pFq({r+1.0, r}, {i+r+1.0}, ((e+d-1.0)/(d-1.0)));
        } catch (...) {
          Rcpp::stop("Error");
        }
      }
    }
  }
  return(hSeq);
}

// not exported
std::vector<mpfr_20> aux_seq_death_n0(mpfr_20 e, mpfr_20 w, mpfr_20 d, int n, int n0=0) {
  std::vector<mpfr_20> hSeq(n+1-n0);
  int i;
  mpfr_20 r=1.0/w;
  if (e==1) {
    if (n0==0) {hSeq[0]=-1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {r+static_cast<mpfr_20>(2.0)}, d)*boost::math::beta(r,2.0);};
    for (i=std::max(n0,1); i<=n; ++i) {
      hSeq[i-n0]=pow(1.0-d,2)*r*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(i+1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+i+1.0)}, d)*boost::math::beta(i,r+1.0);
    };
  } else {
    if (d>=(1.0-2*e)) {
      if (n0==0) {hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, d)+(1.0-e-d)/e*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, ((e+d-1.0)/e)));};
      for (i=std::max(n0,1); i<=n; ++i) {
        hSeq[i-n0] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(i+1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(i+r+1.0)}, ((e+d-1.0)/e));
      }
    } else {
      if (n0==0) {hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {static_cast<mpfr_20>(r+2.0)}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, d)+(1.0-e-d)/e*(e/(1-d))*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(1.0)}, {static_cast<mpfr_20>(r+2.0)}, ((e+d-1.0)/(d-1.0))));};
      for (i=std::max(n0,1); i<=n; ++i) {
        hSeq[i-n0] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*pow(e/(1.0-d),r+1)*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(r+1.0), r}, {static_cast<mpfr_20>(i+r+1.0)}, ((e+d-1.0)/(d-1.0)));
      }
    }
  }
  return(hSeq);
}

// [[Rcpp::export]]
std::vector<double> aux_seq_death_ext_n0(double e, double w, double d, int n, bool boost=false, int n0=0) {
  // NumericVector aux_seq_death_ext(double e, double w, double d, int n, bool boost) {
  std::vector<double> hSeq(n+1-n0);
  int i;
  if (boost==false) {
    std::vector<double> seq(n+1-n0);
    seq=aux_seq_death_n0(e, w, d, n, n0);
    for (i=0;i<=(n-n0);++i) {
      hSeq[i]=seq[i];
    };
    return hSeq;
  } else {
    std::vector<mpfr_20> seq(n+1-n0);
    seq=aux_seq_death_n0(static_cast<mpfr_20>(e), static_cast<mpfr_20>(w), static_cast<mpfr_20>(d), n, n0);
    for (i=0;i<=(n-n0);++i) {
      hSeq[i]=seq[i].convert_to<double>();
    };
    return hSeq;
  }
}

// not exported
void prob_ld(const double &m, int &n0, int &n, std::vector<double> &seq, std::vector<double> &prob){
  int i,j;
  prob.resize(n+1);
  
  if (n0==0) {prob[0]=exp(m*seq[0]);};
  for(i=std::max(n0,1);i<=n;++i) {
    for(j=1;j<=i;++j){
      prob[i]+=j*seq[j]*prob[i-j];
    };
    prob[i]*=m/i;
  };
}

// not exported
void prob_ld(const double &m, int &n0, int &n, std::vector<double> &seq, std::vector<mpfr_20> &prob){
  int i,j;
  prob.resize(n+1);
  
  if (n0==0) {prob[0]=exp(m*seq[0]);};
  for(i=std::max(n0,1);i<=n;++i) {
    for(j=1;j<=i-1;++j){
      prob[i]+=j*seq[j]*prob[i-j];
    };
    prob[i]*=m/i;
  };
}

// not exported
void xi_seq(const double &A, int &n0, int &n, std::vector<double> &seq, std::vector<double> &xi){
  int i,j;
  xi.resize(n+1);
  
  if (n0==0) {xi[0]=log(1.0-A*seq[0]);};
  if (n0==0 || n0==1) xi[1-n0]=-A*seq[1]/(1.0-A*seq[0]);
  for(i=std::max(n0,2);i<=n;++i) {
    for(j=1;j<=i-1;++j) {
      xi[i]+=j*xi[j]*seq[i-j];
    };
    xi[i]=-A*(i*seq[i]-xi[i])/(i*(1.0-A*seq[0]));
  };
}

// not exported
void xi_seq(const mpfr_20 &A, int &n0, int &n, std::vector<double> &seq, std::vector<mpfr_20> &xi){
  int i,j;
  xi.resize(n+1);
  
  if (n0==0) {xi[0]=log(1.0-A*seq[0]);};
  if (n0==0 || n0==1) xi[1-n0]=-A*seq[1]/(1.0-A*seq[0]);
  for(i=std::max(n0,2);i<=n;++i) {
    for(j=1;j<=i-1;++j) {
      xi[i]+=j*xi[j]*seq[i-j];
    };
    xi[i]=-A*(i*seq[i]-xi[i])/(i*(1.0-A*seq[0]));
  };
}

// not exported
void prob_b0(const double &A, const double &k, int &n0, int &n, const double &seq0, std::vector<double> &xi, std::vector<double> &prob){
  int i,j;
  prob.resize(n+1);
  
  if (n0==0) {prob[0]=1.0/pow((1.0-A*seq0),k);};
  for(i=std::max(1,n0);i<=n;++i) {
    for(j=1;j<=i;++j){
      prob[i] += j*xi[j]*prob[i-j];
    };
    prob[i]*=(-k/i);
  };
}

// not exported
void prob_b0(const mpfr_20 &A, const mpfr_20 &k, int &n0, int &n, const double &seq0, std::vector<mpfr_20> &xi, std::vector<mpfr_20> &prob){
  int i,j;
  prob.resize(n+1);
  
  if (n0==0) {prob[0]=1.0/pow((1.0-A*seq0),k);};
  for(i=std::max(1,n0);i<=n;++i) {
    for(j=1;j<=i;++j){
      prob[i] += j*xi[j]*prob[i-j];
    };
    prob[i]*=(-k/i);
  };
}

// not exported
void prob_pois(const double &poisson, int &n0, int &n, std::vector<double> &probp){
  probp.resize(n+1);
  
  if (n0==0) {probp[0]=exp(-poisson);};
  for(int i=std::max(1,n0);i<=n;++i) {
    probp[i]=probp[i-1]*poisson/i;
  };
}

// not exported
void prob_pois(const double &poisson, int &n0, int &n, std::vector<mpfr_20> &probp){
  probp.resize(n+1);
  
  if (n0==0) {probp[0]=exp(-poisson);};
  for(int i=std::max(1,n0);i<=n;++i) {
    probp[i]=probp[i-1]*poisson/i;
  };
}

// not exported
void convolute(int &n0, int &n, std::vector<double> &probc, std::vector<double> &prob, std::vector<double> &probp){
  int i,j;
  probc.resize(n+1);
  
  for(i=n0+1;i<=n+1;++i){
    probc[i-1]=0;
    for(j=1;j<=i;++j){
      probc[i-1]+=prob[j-1]*probp[i-j];
    };
  };
}

// not exported
void convolute(int &n0, int &n, std::vector<mpfr_20> &probc, std::vector<mpfr_20> &prob, std::vector<mpfr_20> &probp){
  int i,j;
  probc.resize(n+1);
  
  for(i=n0+1;i<=n+1;++i){
    probc[i-1]=0;
    for(j=1;j<=i;++j){
      probc[i-1]+=prob[j-1]*probp[i-j];
    };
  };
}

// [[Rcpp::export]]
Rcpp::List prob_mutations_n0_boost(double m, double e=1, double w=1, double cv=0, double death=0, double lag=0, double phi=0, double poisson=0,
                                   double cdf=0.99, int maxiter=10){
  int n=5000;
  int n0=0;
  int iter=1;
  std::vector<double> seq;
  NumericVector seq1(n-n0+1);
  std::vector<mpfr_20> prob(n-n0+1), probp(n-n0+1), probc(n-n0+1), xi(n-n0+1);
  mpfr_20 sumprob=0;
  
  Rcpp::Environment mlemur = Rcpp::Environment::namespace_env("mlemur");
  Rcpp::Function auxseq = mlemur["aux.seq.n0"];
  
  do {
    checkUserInterrupt();
    
    seq1=auxseq(e, w, death, lag, phi, n0, n);
    seq.insert(std::end(seq), std::begin(seq1), std::end(seq1));
    
    if (cv==0) {
      prob_ld(m, n0, n, seq, prob);
    } else {
      mpfr_20 k=1/cv/cv;
      mpfr_20 A=m/k;
      xi_seq(A, n0, n, seq, xi);
      prob_b0(A, k, n0, n, seq[0], xi, prob);
    };
    
    if (poisson!=0) {
      prob_pois(poisson, n0, n, probp);
      convolute(n0, n, probc, prob, probp);
    } else {
      probc.resize(n+1);
      for (int i=n0;i<=n;++i){
        probc[i]=prob[i];
      };
    };
    
    for (int i=n0;i<=n;++i){
      sumprob+=probc[i];
    };

    n0=n+1;
    n+=5000;
    
    iter++;
  } while (sumprob<cdf && iter<maxiter);
  
  NumericVector probres(n+1);
  for (int i=0;i<=n;++i){
    probres[i]=static_cast<double>(probc[i]);
  };
  
  return Rcpp::List::create(Rcpp::Named("seq") = seq,
                            Rcpp::Named("prob") = probres,
                            Rcpp::Named("boost") = true);
}

// [[Rcpp::export]]
Rcpp::List prob_mutations_n0(double m, double e=1, double w=1, double cv=0, double death=0, double lag=0, double phi=0, double poisson=0,
                             double cdf=0.99, int maxiter=10){
  int n=5000;
  int n0=0;
  int iter=1;
  std::vector<double> seq;
  NumericVector seq1(n-n0+1);
  std::vector<double> prob(n-n0+1), probp(n-n0+1), probc(n-n0+1), xi(n-n0+1);
  double sumprob=0;
  
  Rcpp::Environment mlemur = Rcpp::Environment::namespace_env("mlemur");
  Rcpp::Function auxseq = mlemur["aux.seq.n0"];
  
  do {
    checkUserInterrupt();
    
    seq1=auxseq(e, w, death, lag, phi, n0, n);
    seq.insert(std::end(seq), std::begin(seq1), std::end(seq1));
    
    if (cv==0) {
      prob_ld(m, n0, n, seq, prob);
    } else {
      double k=1/cv/cv;
      double A=m/k;
      xi_seq(A, n0, n, seq, xi);
      prob_b0(A, k, n0, n, seq[0], xi, prob);
    };
    
    if (poisson!=0) {
      prob_pois(poisson, n0, n, probp);
      convolute(n0, n, probc, prob, probp);
    } else {
      probc.resize(n+1);
      for (int i=n0;i<=n;++i){
        probc[i]=prob[i];
      };
    };
    
    for (int i=n0;i<=n;++i){
      sumprob+=probc[i];
    };
    if (!std::isfinite(sumprob) || (std::abs(sumprob) < std::numeric_limits<double>::epsilon())) {
      return(prob_mutations_n0_boost(m, e, w, cv, death, lag, phi, poisson, cdf, maxiter));
    };
    
    n0=n+1;
    n+=5000;
    
    iter++;
  } while (sumprob<cdf && iter<maxiter);
  
  return Rcpp::List::create(Rcpp::Named("seq") = seq,
                            Rcpp::Named("prob") = probc,
                            Rcpp::Named("boost") = false);
}

/// HELPER FUNCTIONS FOR ONE-PARAMETER ESTIMATES - SAMPLE SIZE AND POWER EST ///

void derivatives_power(double &U, double &J, double &loglik, const double &m, std::vector<double> &prob0, std::vector<double> &seq, const double &k, const double &poisson, const bool &fisher=true){
  int i;
  int samplesize=prob0.size();
  std::vector<double> prob(samplesize),prob1(samplesize),prob2(samplesize),subscore(samplesize);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, samplesize);
    } else {
      std::vector<double> prob_p(samplesize),prob_lc(samplesize);
      prob_p=prob_pois(poisson, samplesize);
      prob_lc=prob_ld(m, seq, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
    prob1=prob_ld_deriv(seq, prob, samplesize);
    if (fisher) {prob2=prob_ld_deriv(seq, prob1, samplesize);};
  } else {
    std::vector<double> xi(samplesize);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, samplesize);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, samplesize);
      prob1=prob_b0_deriv1(A, k, seq, xi, samplesize);
      if (fisher) {prob2=prob_b0_deriv2(A, k, seq, xi, samplesize);};
    } else {
      std::vector<double> prob_p(samplesize),prob_lc(samplesize),prob_lc_1(samplesize),prob_lc_2(samplesize);
      prob_p=prob_pois(poisson, samplesize);
      prob_lc=prob_b0(A, k, seq0, xi, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
      prob_lc_1=prob_b0_deriv1(A, k, seq, xi, samplesize);
      prob1=prob_ld_deriv(prob_p,prob_lc_1,samplesize);
      if (fisher) {
        prob_lc_2=prob_b0_deriv2(A, k, seq, xi, samplesize);
        prob2=prob_ld_deriv(prob_p,prob_lc_2,samplesize);
      };
    };
  };
  
  U=0;
  J=0;
  loglik=0;
  
  for(i=0;i<=samplesize-1;++i){
    subscore[i]=(prob1[i]/prob[i]);
    U+=subscore[i]*prob0[i];
    if (fisher) {J+=(subscore[i]*subscore[i]-prob2[i]/prob[i])*prob0[i];};
    loglik+=log(prob[i])*prob0[i];
  };
}

void derivatives_power_boost(double &U, double &J, double &loglik, const double &m, std::vector<double> &prob0, std::vector<double> &seq, const double &k, const double &poisson, const bool &fisher=true){
  int i;
  int samplesize=prob0.size();
  std::vector<mpfr_20> prob(samplesize),prob1(samplesize),prob2(samplesize),subscore(samplesize);
  mpfr_20 xm=m,xk=k,xpoisson=poisson,xU=0,xJ=0,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, samplesize);
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_ld(xm, seq, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
    prob1=prob_ld_deriv(seq, prob, samplesize);
    if (fisher) {prob2=prob_ld_deriv(seq, prob1, samplesize);};
  } else {
    std::vector<mpfr_20> xi(samplesize);
    mpfr_20 seq0=seq[0];
    mpfr_20 A=m/k;
    xi=xi_seq(A, seq, samplesize);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, samplesize);
      prob1=prob_b0_deriv1(A, k, seq, xi, samplesize);
      if (fisher) {prob2=prob_b0_deriv2(A, k, seq, xi, samplesize);};
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize),prob_lc_1(samplesize),prob_lc_2(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_b0(A, k, seq0, xi, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
      prob_lc_1=prob_b0_deriv1(A, k, seq, xi, samplesize);
      prob1=prob_ld_deriv(prob_p,prob_lc_1,samplesize);
      if (fisher) {
        prob_lc_2=prob_b0_deriv2(A, k, seq, xi, samplesize);
        prob2=prob_ld_deriv(prob_p,prob_lc_2,samplesize);
      };
    };
  };
  
  xU=0;
  xJ=0;
  xloglik=0;
  
  for(i=0;i<=samplesize-1;++i){
    subscore[i]=(prob1[i]/prob[i]);
    xU+=subscore[i]*prob0[i];
    if (fisher) {xJ+=(subscore[i]*subscore[i]-prob2[i]/prob[i])*prob0[i];};
    xloglik+=log(prob[i])*prob0[i];
  };
  
  U=static_cast<double>(xU);
  J=static_cast<double>(xJ);
  loglik=static_cast<double>(xloglik);
}

void derivatives_power_boost(double &U, double &J, double &loglik, const double &m, std::vector<mpfr_20> &prob0, std::vector<double> &seq, const double &k, const double &poisson, const bool &fisher=true){
  int i;
  int samplesize=prob0.size();
  std::vector<mpfr_20> prob(samplesize),prob1(samplesize),prob2(samplesize),subscore(samplesize);
  mpfr_20 xm=m,xk=k,xpoisson=poisson,xU=0,xJ=0,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, samplesize);
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_ld(xm, seq, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
    prob1=prob_ld_deriv(seq, prob, samplesize);
    if (fisher) {prob2=prob_ld_deriv(seq, prob1, samplesize);};
  } else {
    std::vector<mpfr_20> xi(samplesize);
    mpfr_20 seq0=seq[0];
    mpfr_20 A=m/k;
    xi=xi_seq(A, seq, samplesize);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, samplesize);
      prob1=prob_b0_deriv1(A, k, seq, xi, samplesize);
      if (fisher) {prob2=prob_b0_deriv2(A, k, seq, xi, samplesize);};
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize),prob_lc_1(samplesize),prob_lc_2(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_b0(A, k, seq0, xi, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
      prob_lc_1=prob_b0_deriv1(A, k, seq, xi, samplesize);
      prob1=prob_ld_deriv(prob_p,prob_lc_1,samplesize);
      if (fisher) {
        prob_lc_2=prob_b0_deriv2(A, k, seq, xi, samplesize);
        prob2=prob_ld_deriv(prob_p,prob_lc_2,samplesize);
      };
    };
  };
  
  xU=0;
  xJ=0;
  xloglik=0;
  
  for(i=0;i<=samplesize-1;++i){
    subscore[i]=(prob1[i]/prob[i]);
    xU+=subscore[i]*prob0[i];
    if (fisher) {xJ+=(subscore[i]*subscore[i]-prob2[i]/prob[i])*prob0[i];};
    xloglik+=log(prob[i])*prob0[i];
  };
  
  U=static_cast<double>(xU);
  J=static_cast<double>(xJ);
  loglik=static_cast<double>(xloglik);
}

// not exported
std::vector<double> expected_logliks(double m, double &k, double &poisson, std::vector<double> &seq, std::vector<double> &prob0){
  double E1=0, E2=0;
  int i;
  int samplesize=prob0.size();
  std::vector<double> prob(samplesize);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, samplesize);
    } else {
      std::vector<double> prob_p(samplesize),prob_lc(samplesize);
      prob_p=prob_pois(poisson, samplesize);
      prob_lc=prob_ld(m, seq, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  } else {
    std::vector<double> xi(samplesize);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, samplesize);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, samplesize);
    } else {
      std::vector<double> prob_p(samplesize),prob_lc(samplesize),prob_lc_1(samplesize),prob_lc_2(samplesize);
      prob_p=prob_pois(poisson, samplesize);
      prob_lc=prob_b0(A, k, seq0, xi, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    E1+=log(prob0[i])*prob0[i];
    E2+=log(prob[i])*prob0[i];
  };
  
  return std::vector<double>{E1, E2};
}

// not exported
std::vector<double> expected_logliks(double m, double &k, double &poisson, std::vector<double> &seq, std::vector<mpfr_20> &prob0){
  mpfr_20 E1=0, E2=0;
  int i;
  int samplesize=prob0.size();
  std::vector<double> prob(samplesize);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, samplesize);
    } else {
      std::vector<double> prob_p(samplesize),prob_lc(samplesize);
      prob_p=prob_pois(poisson, samplesize);
      prob_lc=prob_ld(m, seq, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  } else {
    std::vector<double> xi(samplesize);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, samplesize);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, samplesize);
    } else {
      std::vector<double> prob_p(samplesize),prob_lc(samplesize),prob_lc_1(samplesize),prob_lc_2(samplesize);
      prob_p=prob_pois(poisson, samplesize);
      prob_lc=prob_b0(A, k, seq0, xi, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    E1+=log(prob0[i])*prob0[i];
    E2+=log(prob[i])*prob0[i];
  };
  
  return std::vector<double>{static_cast<double>(E1), static_cast<double>(E2)};
}

// not exported
std::vector<double> expected_logliks_boost(double m, double &k, double &poisson, std::vector<double> &seq, std::vector<double> &prob0){
  mpfr_20 E1=0, E2=0;
  int i;
  int samplesize=prob0.size();
  std::vector<mpfr_20> prob(samplesize);
  mpfr_20 xm=m,xk=k,xpoisson=poisson,xU=0,xJ=0,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, samplesize);
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_ld(xm, seq, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  } else {
    std::vector<mpfr_20> xi(samplesize);
    mpfr_20 seq0=seq[0];
    mpfr_20 A=m/k;
    xi=xi_seq(A, seq, samplesize);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, samplesize);
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize),prob_lc_1(samplesize),prob_lc_2(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_b0(A, k, seq0, xi, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    E1+=log(prob0[i])*prob0[i];
    E2+=log(prob[i])*prob0[i];
  };
  
  return std::vector<double>{static_cast<double>(E1), static_cast<double>(E2)};
}

// not exported
std::vector<double> expected_logliks_boost(double m, double &k, double &poisson, std::vector<double> &seq, std::vector<mpfr_20> &prob0){
  mpfr_20 E1=0, E2=0;
  int i;
  int samplesize=prob0.size();
  std::vector<mpfr_20> prob(samplesize);
  mpfr_20 xm=m,xk=k,xpoisson=poisson,xU=0,xJ=0,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, samplesize);
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_ld(xm, seq, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  } else {
    std::vector<mpfr_20> xi(samplesize);
    mpfr_20 seq0=seq[0];
    mpfr_20 A=m/k;
    xi=xi_seq(A, seq, samplesize);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, samplesize);
    } else {
      std::vector<mpfr_20> prob_p(samplesize),prob_lc(samplesize),prob_lc_1(samplesize),prob_lc_2(samplesize);
      prob_p=prob_pois(xpoisson, samplesize);
      prob_lc=prob_b0(A, k, seq0, xi, samplesize);
      prob=prob_ld_deriv(prob_p,prob_lc,samplesize);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    E1+=log(prob0[i])*prob0[i];
    E2+=log(prob[i])*prob0[i];
  };
  
  return std::vector<double>{static_cast<double>(E1), static_cast<double>(E2)};
}

// [[Rcpp::export]]
std::vector<double> optim_m_from_probs(double &current_m, double &lower_m, double &upper_m, double &R,
                                       double &k1, double &k2, double &poisson1, double &poisson2,
                                       std::vector<double> &seq1, std::vector<double> &seq2,
                                       std::vector<double> &prob1, std::vector<double> &prob2,
                                       double &m1, double &m2, bool &prob1_boost, bool &prob2_boost,
                                       bool verbose=false){
  int iter=0;
  double new_m,new_U,new_J,new_loglik,current_U,current_J,current_loglik,lower_U,lower_J,lower_loglik,upper_U,upper_J,upper_loglik,
  U1,J1,loglik1,U2,J2,loglik2;
  
  std::vector<mpfr_20> prob1x(seq1.size()), prob2x(seq2.size());
  
  double cv1, cv2;
  if (k1 == 0) {cv1 = 0;} else {cv1 = 1/sqrt(k1);}
  if (k2 == 0) {cv2 = 0;} else {cv2 = 1/sqrt(k2);}
  
  if (prob1_boost) {prob1x = calc_probs_0(static_cast<mpfr_20>(m1), cv1, seq1, seq1.size(), poisson1);};
  if (prob2_boost) {prob1x = calc_probs_0(static_cast<mpfr_20>(m2), cv2, seq2, seq2.size(), poisson2);};
  
  if (!prob1_boost) {derivatives_power(U1, J1, loglik1, upper_m, prob1, seq1, k1, poisson1, true);} else {derivatives_power_boost(U1, J1, loglik1, upper_m, prob1x, seq1, k1, poisson1, true);}
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_power_boost(U1, J1, loglik1, upper_m, prob1, seq1, k1, poisson1, true);};
  if (!prob2_boost) {derivatives_power(U2, J2, loglik2, R*upper_m, prob2, seq2, k2, poisson2, true);} else {derivatives_power_boost(U2, J2, loglik2, R*upper_m, prob2x, seq2, k2, poisson2, true);}
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_power_boost(U2, J2, loglik2, R*upper_m, prob2, seq2, k2, poisson2, true);};
  upper_U=U1+R*U2;
  upper_J=J1+R*R*J2;
  upper_loglik=loglik1+loglik2;
  if (std::abs(upper_U)<1e-6) {return std::vector<double>{upper_m,upper_loglik};};
  
  if (!prob1_boost) {derivatives_power(U1, J1, loglik1, lower_m, prob1, seq1, k1, poisson1, true);} else {derivatives_power_boost(U1, J1, loglik1, lower_m, prob1x, seq1, k1, poisson1, true);}
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_power_boost(U1, J1, loglik1, lower_m, prob1, seq1, k1, poisson1, true);};
  if (!prob2_boost) {derivatives_power(U2, J2, loglik2, R*lower_m, prob2, seq2, k2, poisson2, true);} else {derivatives_power_boost(U2, J2, loglik2, R*lower_m, prob2x, seq2, k2, poisson2, true);}
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_power_boost(U2, J2, loglik2, R*lower_m, prob2, seq2, k2, poisson2, true);};
  lower_U=U1+R*U2;
  lower_J=J1+R*R*J2;
  lower_loglik=loglik1+loglik2;
  if (std::abs(lower_U)<1e-6) {return std::vector<double>{lower_m,lower_loglik};};
  
  while(upper_U>0){
    checkUserInterrupt();
    lower_m=upper_m;lower_U=upper_U;lower_J=upper_J;lower_loglik=upper_loglik;
    upper_m*=5;
    if (!prob1_boost) {derivatives_power(U1, J1, loglik1, upper_m, prob1, seq1, k1, poisson1, true);} else {derivatives_power_boost(U1, J1, loglik1, upper_m, prob1x, seq1, k1, poisson1, true);}
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_power_boost(U1, J1, loglik1, upper_m, prob1, seq1, k1, poisson1, true);};
    if (!prob2_boost) {derivatives_power(U2, J2, loglik2, R*upper_m, prob2, seq2, k2, poisson2, true);} else {derivatives_power_boost(U2, J2, loglik2, R*upper_m, prob2x, seq2, k2, poisson2, true);}
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_power_boost(U2, J2, loglik2, R*upper_m, prob2, seq2, k2, poisson2, true);};
    upper_U=U1+R*U2;
    upper_J=J1+R*R*J2;
    upper_loglik=loglik1+loglik2;
    if (std::abs(upper_U)<1e-6) {return std::vector<double>{upper_m,upper_loglik};};
  };
  while(lower_U<0){
    checkUserInterrupt();
    upper_m=lower_m;upper_U=lower_U;upper_J=lower_J;upper_loglik=lower_loglik;
    lower_m/=10;
    if (lower_m<1e-24) {return {0};};
    if (!prob1_boost) {derivatives_power(U1, J1, loglik1, lower_m, prob1, seq1, k1, poisson1, true);} else {derivatives_power_boost(U1, J1, loglik1, lower_m, prob1x, seq1, k1, poisson1, true);}
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_power_boost(U1, J1, loglik1, lower_m, prob1, seq1, k1, poisson1, true);};
    if (!prob2_boost) {derivatives_power(U2, J2, loglik2, R*lower_m, prob2, seq2, k2, poisson2, true);} else {derivatives_power_boost(U2, J2, loglik2, R*lower_m, prob2x, seq2, k2, poisson2, true);}
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_power_boost(U2, J2, loglik2, R*lower_m, prob2, seq2, k2, poisson2, true);};
    lower_U=U1+R*U2;
    lower_J=J1+R*R*J2;
    lower_loglik=loglik1+loglik2;
    if (std::abs(lower_U)<1e-6) {return std::vector<double>{lower_m,lower_loglik};};
  };
  
  if(verbose) {Rcout << "boundaries:: m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "boundaries:: U: " << lower_U << " " << upper_U << "\n";};
  if(verbose) {Rcout << "boundaries:: loglik: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m<lower_m) || (current_m>upper_m)) {current_m=(lower_m+upper_m)/2;};
  
  if (!prob1_boost) {derivatives_power(U1, J1, loglik1, current_m, prob1, seq1, k1, poisson1, true);} else {derivatives_power_boost(U1, J1, loglik1, current_m, prob1x, seq1, k1, poisson1, true);}
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_power_boost(U1, J1, loglik1, current_m, prob1, seq1, k1, poisson1, true);};
  if (!prob2_boost) {derivatives_power(U2, J2, loglik2, R*current_m, prob2, seq2, k2, poisson2, true);} else {derivatives_power_boost(U2, J2, loglik2, R*current_m, prob2x, seq2, k2, poisson2, true);}
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_power_boost(U2, J2, loglik2, R*current_m, prob2, seq2, k2, poisson2, true);};
  current_U=U1+R*U2;
  current_J=J1+R*R*J2;
  current_loglik=loglik1+loglik2;
  
  if(verbose) {Rcout<< "U: " << current_U << " J: " << current_J << " loglik: " << current_loglik << " m: " << current_m << "\n\n";};
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {return {-1.0};};
  
  while(std::abs(current_U)>1.0e-6){
    checkUserInterrupt();
    iter++;
    if (iter>50) {return {-1.0};};
    
    new_m=current_m+current_U/current_J;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    };
    
    if (!prob1_boost) {derivatives_power(U1, J1, loglik1, new_m, prob1, seq1, k1, poisson1, true);} else {derivatives_power_boost(U1, J1, loglik1, new_m, prob1x, seq1, k1, poisson1, true);}
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_power_boost(U1, J1, loglik1, new_m, prob1, seq1, k1, poisson1, true);};
    if (!prob2_boost) {derivatives_power(U2, J2, loglik2, R*new_m, prob2, seq2, k2, poisson2, true);} else {derivatives_power_boost(U2, J2, loglik2, R*new_m, prob2x, seq2, k2, poisson2, true);}
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_power_boost(U2, J2, loglik2, R*new_m, prob2, seq2, k2, poisson2, true);};
    new_U=U1+R*U2;
    new_J=J1+R*R*J2;
    new_loglik=loglik1+loglik2;
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return {-1.0};};
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if (new_U<0){
      upper_m=new_m;upper_U=new_U;upper_J=new_J;upper_loglik=new_loglik;
    } else if (new_U>0){
      lower_m=new_m;lower_U=new_U;lower_J=new_J;lower_loglik=new_loglik;
    };
    current_m=new_m;current_U=new_U;current_J=new_J;current_loglik=new_loglik;
  };
  
  std::vector<double> ret(6), loglika(2), loglikb(2);
  if (!prob1_boost) {
    loglika = expected_logliks(current_m, k1, poisson1, seq1, prob1);
    if (std::any_of(loglika.begin(), loglika.end(), [](double i){return(!std::isfinite(i));})) {loglika = expected_logliks_boost(current_m, k1, poisson1, seq1, prob1);}
  } else {
    loglika = expected_logliks(current_m, k1, poisson1, seq1, prob1x);
    if (std::any_of(loglika.begin(), loglika.end(), [](double i){return(!std::isfinite(i));})) {loglika = expected_logliks_boost(current_m, k1, poisson1, seq1, prob1x);}
  };
  if (!prob2_boost){
    loglikb = expected_logliks(current_m*R, k2, poisson2, seq2, prob2);
    if (std::any_of(loglikb.begin(), loglikb.end(), [](double i){return(!std::isfinite(i));})) {loglikb = expected_logliks_boost(R*current_m, k2, poisson2, seq2, prob2);}
  } else {
    loglikb = expected_logliks(current_m*R, k2, poisson2, seq2, prob2x);
    if (std::any_of(loglikb.begin(), loglikb.end(), [](double i){return(!std::isfinite(i));})) {loglikb = expected_logliks_boost(R*current_m, k2, poisson2, seq2, prob2x);}
  };
  ret = {current_m, current_loglik, loglika[0], loglika[1], loglikb[0], loglikb[1]};
  return ret;
}


