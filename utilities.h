//
// Created by Sindre Bakke Øyen on 14.02.2018.
//

#ifndef MASTERPROJECTCPP_READUTILS_H
#define MASTERPROJECTCPP_READUTILS_H
/* Built-in header files */
#include <stdio.h>                      /* Used for functions such as strsep and strtod */
#include <string.h>                     /* Used for strings */
#include <stdlib.h>                     /* Used for functions such as printf */

/* External library header files */
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>     /* access to serial N_Vector  */
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix  */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* User-defined header files */

/* DECLARATIONS */
void getRowsAndCols(const char *filename, size_t *rows, size_t *cols);
void getDistributions(const char *filename, gsl_matrix *fv, gsl_vector *r, gsl_vector *t,
                      size_t rows, size_t cols);
void rescaleInitial(gsl_matrix *fv, gsl_vector *r, realtype phi);

#endif //MASTERPROJECT_READUTILS_H
