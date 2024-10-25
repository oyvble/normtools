//#include <cmath> //includes special functions
#include <RcppArmadillo.h>
using namespace arma;

extern "C" {
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
 
void EM_estLoop(double *xptr, double *lambda, double *mu, double *sigmaSq, double *p1, int *n, double *epsilon) {
  //Rcpp::NumericVector retvec(4); //init return vector
  static const double pi = 3.14159265;
  
  //copy variables
  vec *x = new vec(xptr, *n, false);
  vec x2 = (*x)%(*x);  //Calc elementwise/Schur product and keeping long vector
  vec z; 
  long double sumX = sum(*x);
  long double con,mu2,sumZ,sumXZ,sum1Z; //help variables
  long double lambda_new,mu_new,sigmaSq_new,p1_new; //declar variables
  double log2pi = log(2*pi); //obtain constant
  bool diff = true; //whether there are any difference
  while (diff) { 
    mu2 = (*mu)*(*mu); //mu^2
    con = 0.5*(mu2/ *sigmaSq + log(*sigmaSq) + log2pi) + log(*lambda); //constant to update
    z = -0.5*x2/ *sigmaSq + (*x)*( *mu / *sigmaSq + *lambda) - con; //calulcate logged expression
    z = *p1/(*p1+(1-*p1)*exp(z));
    
    sumZ = sum(z);
    sumXZ = sum((*x)%z);
    sum1Z = (*n - sumZ); //sum(1-z);
    
    p1_new = sumZ/ *n;
    lambda_new = sumZ/sumXZ;
    mu_new = (sumX - sumXZ)/sum1Z; //sum((1-z)*x)/sum(1-z);
    
    vec tmp = *x-mu_new;
    sigmaSq_new = sum( tmp%tmp%(1-z) )/sum1Z;
    
    //diff <- !(all(c(abs(lambda_new-*lambda),abs(mu_new-*mu),abs(sigmaSq_new-*sigmaSq),abs(p1_new-*p1))<epsilon))
    diff = (abs(lambda_new-*lambda)>= epsilon[0]) || (abs(mu_new-*mu)>=epsilon[1]) || (abs(sigmaSq_new-*sigmaSq)>=epsilon[2]) || (abs(p1_new-*p1)>=epsilon[3]); //check if any above epsilon in difference
    
    //Update variables:
    *lambda = lambda_new; 
    *mu = mu_new; 
    *sigmaSq = sigmaSq_new; 
    *p1 = p1_new;
  }
   delete x;
}

} //end external
