#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <armadillo>




arma::cx_dvec expmv(std::complex<double> t, arma::sp_cx_dmat& A, arma::cx_dvec& v, double tol = 1.e-7, unsigned int m = 30);

bool rsvd(arma::cx_dmat& U, arma::dvec& s, arma::cx_dmat& V, arma::cx_dmat& D, unsigned int k);



