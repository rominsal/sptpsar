#include <RcppEigen.h>
// using namespace Eigen;
typedef Eigen::Map<Eigen::MatrixXd>  MapMatd;
typedef Eigen::Map<Eigen::MatrixXi>  MapMati;
typedef Eigen::Map<Eigen::VectorXd>  MapVecd;


// [[Rcpp::depends(RcppEigen)]]

// Functions for dense matrices

// Cross-Products

// [[Rcpp::export]]
    Eigen::MatrixXd AtA(const Eigen::MatrixXd & x) {
    int    n(x.cols());
    Eigen::MatrixXd xtx = Eigen::MatrixXd(n,n).setZero().
    selfadjointView<Eigen::Lower>().rankUpdate(x.adjoint());
  // Eigen::MatrixXd m = x.transpose() * x;
    return xtx;
  }

// [[Rcpp::export]]
    Eigen::MatrixXd AAt(const Eigen::MatrixXd & x) {
              int    n(x.cols());
              Eigen::MatrixXd xxt = Eigen::MatrixXd(n,n).setZero().
                         selfadjointView<Eigen::Lower>().rankUpdate(x);
      //Eigen::MatrixXd m = x * x.transpose();
      return xxt;
    }

// [[Rcpp::export]]
    Eigen::MatrixXd lltA(const Eigen::MatrixXd & x) {
      const Eigen::LLT<Eigen::MatrixXd> llt(x);
      Eigen::MatrixXd R = llt.matrixU();
      return R;
    }

// Inverse for Positive Definite Matrix (Cholesky)
// [[Rcpp::export]]
    Eigen::MatrixXd lltAinv(const Eigen::MatrixXd & x) {
      const int     p(x.cols());
      const Eigen::LLT<Eigen::MatrixXd> llt(x);
      Eigen::MatrixXd Ainv = llt.solve(Eigen::MatrixXd::Identity(p, p));
      return Ainv;
  }


// Det and Log-Det for Definite Positive Matrix
// [[Rcpp::export]]
    Rcpp::List detlltA(const Eigen::MatrixXd & x){
                const Eigen::LDLT<Eigen::MatrixXd> ldlt(x);
                const Eigen::VectorXd Dvec(ldlt.vectorD());
                return Rcpp::List::create(
                Rcpp::Named("det") = Dvec.prod(),
                Rcpp::Named("ldet") = Dvec.array().log().sum());
    }


// Inverse for Definite Positive Matrix
// [[Rcpp::export]]
    Eigen::MatrixXd Ainv(const Eigen::MatrixXd & x){
              Eigen::MatrixXd xinv = x.inverse();
              return xinv;
    }


// Functions for sparse matrices

// [[Rcpp::export]]
    Eigen::SparseMatrix<double> AtAsp(const Eigen::SparseMatrix<double>  & x){
      Eigen::SparseMatrix<double> xtx = x.transpose() * x;
      return xtx;
    }

// [[Rcpp::export]]
    Eigen::SparseMatrix<double> lltAsp(const Eigen::SparseMatrix<double> & x) {
      const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> >llt(x);
      Eigen::SparseMatrix<double>  R = llt.matrixU();
      return R;
    }

// Det and Log-Det for Definite Positive Sparse Matrix (Cholesky)
// [[Rcpp::export]]
    Rcpp::List detlltAsp(const Eigen::SparseMatrix<double> & x) {
      const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldlt(x);
      const Eigen::VectorXd Dvec(ldlt.vectorD());
      return Rcpp::List::create(
        Rcpp::Named("det") = Dvec.prod(),
        Rcpp::Named("ldet") = Dvec.array().log().sum());
    }

// Log-Det for Sparse Matrix (LU)
// [[Rcpp::export]]
    double ldetluAsp(const Eigen::SparseMatrix<double> & x) {
      const Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(x) ;
      // const Eigen::VectorXd Dvec(ldlt.vectorD());
      return solver.logAbsDeterminant();
    }
