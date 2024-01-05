#include "HermitianOperator.h"

HermitianOperator::HermitianOperator():arma::sp_cx_fmat(){};

HermitianOperator::HermitianOperator(int basisSize):arma::sp_cx_fmat(basisSize, basisSize){}; 	 

HermitianOperator::HermitianOperator(arma::umat locations, arma::cx_frowvec values, int basisSize):arma::sp_cx_fmat(true,  locations, values, basisSize, basisSize){};

HermitianOperator::HermitianOperator(arma::sp_cx_fmat mat):arma::sp_cx_fmat(mat){};

float HermitianOperator::evaluate(arma::cx_fvec psi)
{
    return static_cast<float>(real(cdot(psi, (arma::sp_cx_fmat)*this * psi)));
}

 arma::cx_fvec HermitianOperator::expmv(arma::cx_float alpha , arma::cx_fvec psi)
 {
    return (psi + alpha * (arma::sp_cx_fmat)*this *psi);
 }