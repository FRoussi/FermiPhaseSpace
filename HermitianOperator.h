#pragma once

#include <armadillo>

class HermitianOperator : public arma::sp_cx_fmat 	
{
public:
    HermitianOperator();	 	 
    HermitianOperator(int basisSize);
    HermitianOperator(arma::umat locations, arma::cx_frowvec values, int basisSize);
    HermitianOperator(arma::sp_cx_fmat mat);

    float evaluate(arma::cx_fvec psi);
    arma::cx_fvec expmv(arma::cx_float alpha , arma::cx_fvec psi);
    
};
