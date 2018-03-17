#include <iostream>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>

/* User-defined header files */
#include "Kernels.h"
#include "utilities.h"

int main() {
    //////////////////////////// DESCRIPTION OF PROGRAM ////////////////////////////
    /* This is the main program that solves the Population Balance Equation (PBE)
     * The SUNMatrices used are stored column-wise, i.e.
     *         time ------------------------------>
     *        -                                   -
     *     x |a(xi0, t0) a(xi0, t1) ... a(xi0, tN) |
     *     i |a(xi1, t0) a(xi1, t1) ... a(xi1, tN) |
     *     | |     .      .                  .     |
     * A = | |     .             .           .     | in R^[MxN]
     *     | |     .                   .     .     |
     *     ↓ |a(xiM, t0) a(xiM, t1) ... a(xiM, tN) |
     *        -                                   -
     * This means we have to retrieve the columns of A
     *      realtype **cols = SM_COLS_D(A);
     * and index it as such:
     *      cols[j][i], where 0 < j < N and 0 < i < M
     * Thus i points to the row, and j points to the column in the matrix A.
     */
    //////////////////////////// DECLARE VARIABLES ////////////////////////////
    /* Set number of non-relevant rows and columns */
    const size_t trashRows = 2;
    const size_t trashCols = 9;

    /* Declare needed variables */
    gsl_matrix *fv;     /* Droplet distributions */
    gsl_vector *r, *t;  /* Size classes and time of measurements */

    size_t rows, cols;  /* Rows and columns of CSV sheet */
    size_t M, N;        /* Rows and columns to preallocate */
    size_t Np = 10;     /* Number of grid points */

    //////////////////////////// EXTRACT FROM CSV ////////////////////////////
    /* Get rows and columns and allocate vector and matrix */
    getRowsAndCols("crudeB.csv", &rows, &cols);
    M = rows - trashRows;           /* Number of rows we will use           */
    N = cols - trashCols;           /* Number of columns we will use        */
    r = gsl_vector_alloc(N);
    t = gsl_vector_alloc(M);
    fv = gsl_matrix_alloc(M, N);

    getDistributions("crudeB.csv", fv, r, t, rows, cols);
    rescaleInitial(fv, r, 0.7e-2);


    Grid g = Grid(Np, 0, 1, 0, 0, 2);
    Fluid disp = Fluid(0.837e3, 22.0e-3, 16.88e-3);
    Fluid cont = Fluid(1.0e3, 1, 1);
    SystemProperties s = SystemProperties(500.0e-6, 725.0e-6, 0.366, 2965, disp);
    Constants c = Constants(1.8190e-12,1.819e-9,1,1.0e3,s, cont, disp);
    Kernels k = Kernels(c,g,s);
    std::cout << k << std::endl;

    // TODO: Create the model to solve
    return 0;
}