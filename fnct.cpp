#include "fnct.h"






/// @brief Copy of the EXPOKIT expmv function
/// @param t double
/// @param A sparse matrix
/// @param v vector
/// @param tol
/// @param m1
/// @return exp(t.A)v
arma::cx_dvec expmv(std::complex<double> t, arma::sp_cx_dmat& A, arma::cx_dvec& v, double tol, unsigned int m1 )
{
    unsigned int n = A.n_rows;
    unsigned int m = std::min(n, m1);
    double anorm = arma::norm(A, "inf");
    unsigned int mxrej = 10;
    double btol = 1.E-7;
    double gamma = 0.9;
    double delta = 1.2;
    unsigned int mb = m;
    unsigned int mx = mb;
    double t_out = std::abs(t);
    int nstep = 0;
    double t_new = 0;
    double t_now = 0;
    double t_step = 0;
    double s = 0;
    double s_error = 0;
    double rndoff = anorm * 1.E-15;

    unsigned int k1 = 2;
    double xm = 1/(double)m;
    double normv = arma::norm(v,2);
    double avnorm = 0;
    double phi1 = 0;
    double phi2 = 0;
    double err_loc = 0;
    double beta = normv;
    double fact = std::pow((m+1)/std::exp(1), m+1) * std::sqrt(2*M_PI*(m+1));
    t_new = (1/anorm) * pow((fact * tol) / 4*beta*anorm, xm);     
    std::complex<double> sgn = t/std::abs(t);
    nstep = 0;

    arma::cx_dvec w = v;
    double hump = normv;

    while(t_now < t_out)
    {
        nstep +=1;
        t_step = std::min(t_out - t_now, t_new);
        arma::cx_dmat V = arma::cx_dmat(n, m+1); 
        arma::cx_dmat H = arma::cx_dmat(m+2, m+2);
        arma::cx_dmat F;

        V.col(0) = (1/beta) * w;
        for (unsigned int j = 0; j<m; j++)
        {
            arma::cx_dvec p = A * V.col(j);
            for (unsigned int i = 0; i < j ; i++)
            {
                H(i,j) = arma::cdot(V.col(i), p);
                p = p - H(i,j)*V.col(i);
            }
            s = arma::norm(p,2);
            if (s<btol)
            {
                k1 = 0;
                mb = j;
                t_step = t_out - t_now;
                break;
            }
            H(j+1, j) = s;
            V.col(j+1) = (1/s) * p;
        }

        if (k1 != 0)
        {
            H(m+1, m) = 1;
            avnorm = arma::norm(A*V.col(m), 2);
        }
        unsigned int ireject = 0;
        while(ireject < mxrej)
        {
            mx = mb + k1;
            F = arma::expmat(sgn * t_step * H.submat(0, 0, mx-1, mx-1));
            if (k1 ==0)
            {
                err_loc = btol;
                break;
            }
            else
            {
                phi1 = std::abs(beta * F(m, 0));
                phi2 = std::abs(beta * F(m+1,0) * avnorm);
                if (phi1 > 10*phi2)
                {
                    err_loc = phi2;
                    xm = 1/(double)m;
                }
                else if (phi1 > phi2)
                {
                    err_loc = (phi1*phi2)/(phi1-phi2);
                    xm = 1/(double)m;
                }
                else
                {
                    err_loc = phi1;
                    xm = 1/(double)(m-1);
                }
            }
            if (err_loc <= delta * t_step * tol)
            {
                break;
            }
            else
            {
                t_step = gamma * t_step * pow((t_step * tol / err_loc), xm);
                s = pow(10, std::floor(std::log10(t_step))-1);
                t_step = std::ceil(t_step/s) * s;
                if (ireject == mxrej)
                {
                    printf("The requested tolerance is too high.");
                    // !!!! THROW AN ERROR HERE
                }
                ireject+=1;
            }
        }
        mx = mb + std::max((unsigned int)0, k1-1);
        w = V.cols(0, mx-1) * (beta * (F.col(0)).subvec(0,mx-1));
        beta = arma::norm(w, 2);
        hump = std::max(hump, beta);

        t_now = t_now + t_step;
        t_new = gamma * t_step * pow((t_step*tol)/err_loc, xm);
        s = pow(10, std::floor(std::log10(t_new))-1);
        t_new = std::ceil(t_new / s) * s;

        err_loc = std::max(err_loc, rndoff);
        s_error = s_error + err_loc;
    }

    return w;

}




/// @brief Copy of the MATLAB rsvd function.  random SVD
/// Extremely fast computation of the truncated Singular Value Decomposition, using
/// randomized algorithms as described in Halko et al. 'finding structure with randomness'.
/// Antoine Liutkus  (c) Inria 2014
/// @param U 
/// @param S matrix of singular values
/// @param V
/// @param D matrix whose svd we want
/// @param K number of compponents to keep
/// @return false if the decomposition failed
bool rsvd(arma::cx_dmat& U, arma::dvec& s, arma::cx_dmat& V, arma::cx_dmat& D, unsigned int K)
{
    bool successOrth = true;
    bool successSVD = true; 
    
                                                    // MATLAB version
    unsigned int M = D.n_rows;
    unsigned int N = D.n_cols;                      //<- [M,N] = arma::size(D);

    unsigned int P = std::min(2*K,N);               //<- P = min(2*K,N);


    arma::dmat X = arma::dmat(N, P, arma::fill::randn); //<- X = randn(N,P);
    arma::cx_dmat Y = D*X;                          //<- Y = A*X;

    arma::cx_dmat W1;
    successOrth = arma::orth(W1, Y);               //<- W1 = orth(Y);

    arma::cx_dmat B = W1.t() * D;                   //<- B = W1'*A;
    
    arma::cx_dmat W2; 

    successSVD = arma::svd_econ(W2, s, V, B);       //<- [W2,S,V] = svd(B,'econ');

    U = W1 * W2;                                    //<- U = W1*W2;

    unsigned int uCols = U.n_cols;

    unsigned int K2 = std::min(K,uCols);            //<- K = min(K, size(U,2))

    U = U.head_cols(K2);                            //<- U = U(:,1:K);
    s =s.head(K2);                                  //<- S = S(1:K,1:K);
    V = V.head_cols(K2);                            //<- V=V(:,1:K);

    return successOrth && successSVD;
} 