//#include <stdio>
#include <iostream>
#include <string>

//#define ARMA_DONT_USE_OPENMP
#include <armadillo>


//#include <filesystem>


#define DEF_MAX_THREAD_NUM 4

int main(int argc, char* argv[])
{



   std::cout << "--- We're here for tests --- " << std::endl;


   arma::arma_version ver;
   std::cout << "== ARMA version: "<< ver.as_string() << std::endl;







// Threads used by openblas :
   //openblas_set_num_threads(1);


   std::cout << "== Double precision :   2^ -" << std::numeric_limits<double>::digits10 << std::endl;


uint n = 3; 

   arma::sp_dmat A =  arma::sprandn(n, n, 0.1);

   arma::cx_dvec x = arma::cx_dvec(3,arma::fill::randn);

   arma::cx_dvec y = A*x;

A.print();

   y.print();
   y.t().print();

   arma::dmat B = (A + A.t()).as_dense();

   B.print();


std::cout << "<yy>" << arma::cdot(y, y) << std::endl;
   std::cout << "<yBy>" << arma::cdot(y, B*y) << std::endl;

   //std::string work = "2Sites.json";
 

   std::cout <<"== Program closing." << std::endl;
   return 0;
}
