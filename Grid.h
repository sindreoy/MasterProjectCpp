//
// Created by Sindre Bakke Øyen on 06.03.2018.
//

#ifndef MASTERPROJECTCPP_GRID_H
#define MASTERPROJECTCPP_GRID_H
/* Built-in header files */
#include <ostream>                      /* Print to console                     */
#include <iomanip>                      /* Manipulate output format             */

/* External library header files */
#include <sundials/sundials_types.h>    /* Datatypes from sundials              */
#include <sundials/sundials_math.h>     /* Math functions, power etc            */
#include <gsl/gsl_vector_double.h>      /* Vectors                              */
#include <gsl/gsl_matrix_double.h>      /* Matrices                             */
#include <gsl/gsl_linalg.h>             /* Linear algebra                       */
#include <gsl/gsl_blas.h>               /* Basic linear algebraic subprograms   */
#include <gsl/gsl_eigen.h>              /* Eigenvectors and values              */

/* User-defined header files */

class Grid {
private:
    gsl_vector *xi, *w;
    gsl_matrix *D, *xipBB, *xipBC, *xippBC;
    size_t N;
    realtype x0, x1, alpha, beta, mu0;
    /* Variables:
     * xi   :: Quadrature points
     * w    :: Quadrature weights
     * D    :: Lagrange derivative matrix
     * xip's:: Interpolated quadrature points for birth terms
     * N    :: Number of grid points
     * x0   :: Left boundary
     * x1   :: Right boundary
     * alpha:: Chooses Jacobi polynomial
     * beta :: Chooses Jacobi polynomial
     * mu0  :: The integral of the weight function in the domain [-1,1] */
public:
    Grid();
    Grid(const Grid &g);
    Grid(size_t N, realtype x0, realtype x1, realtype alpha, realtype beta, realtype mu0);

    void coefs(size_t j, realtype *r);
    void setQuadratureRule();   /* Gauss Lobatto rule */
    void remapGrid();           /* Remaps grid to [x0, x1] domain */
    void setLagrangeDerivativeMatrix();
    void setInterpolatedXis();

    size_t getN() const;
    realtype getX0() const;
    realtype getX1() const;
    realtype getAlpha() const;
    realtype getBeta() const;
    realtype getMu0() const;

    gsl_vector *getXi() const;
    gsl_vector *getW() const;
    gsl_matrix *getD() const;
    gsl_matrix *getXipBB() const;
    gsl_matrix *getXipBC() const;
    gsl_matrix *getXippBC() const;

    friend std::ostream &operator<<(std::ostream &os, const Grid &grid);

    ~Grid();
};

#endif //MASTERPROJECTCPP_GRID_H
