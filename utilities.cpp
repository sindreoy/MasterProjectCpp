//
// Created by Sindre Bakke Øyen on 14.02.2018.
//

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include "utilities.h"

void getRowsAndCols(const char *filename, size_t *rows, size_t *cols){
    /* Input args:
     *  filename    :: pointer to CSV filename to read from
     *  rows        :: total number of rows to count in CSV
     *  cols        :: total number of columns to count in CSV
     * */
    FILE *f = fopen(filename, "r");
    if (f != nullptr) {
        // Initialize variables
        *rows = 0;
        *cols = 0;
        size_t tmp = *cols;
        const char s[2] = ",";
        const int bufSize = 1000;
        char line[bufSize], *tofree;

        // Now loop over lines
        while(fgets(line, sizeof line, f) != nullptr) {
            (*rows)++;
            tmp = 0;
            tofree = strdup(line);
            while ((strsep(&tofree, s)) != nullptr) {
                (tmp)++;
            }
            if (tmp > *cols) {
                *cols = tmp;
            }
        }
        fclose(f);
    } else {
        perror(filename);
    }
}

void getDistributions(const char *filename, gsl_matrix *fv, gsl_vector *r, gsl_vector *t,
                      const size_t rows, const size_t cols){
    /* Input args:
     *  filename    :: pointer to a CSV filename to read from
     *  fv          :: matrix of MxN rows and columns to store measured distributions in
     *  r           :: vector of N columns to store measured droplet size classes in
     *  t           :: vector of M rows to store time of measurement in
     *  rows        :: total number of rows in CSV
     *  cols        :: total number of columns in CSV
     * */
    /* Set number of non-relevant rows and columns */
    const size_t trashRows = 1;
    const size_t trashCols = 9;

    /* Declare needed variables */
    char *hours, *minutes;
    realtype h, m;
    size_t i = 0; // Index variable
    size_t k = 0; // Index variable

    const char s[2] = ",";
    const size_t bufSize = 10*(cols);
    char line[bufSize], *token, *tofree, *time;

    /* Open file and start reading */
    FILE *f = fopen(filename, "r");
    if (f != nullptr) {
        while (fgets(line, sizeof line, f) != nullptr) {
            if (i == 0) {
                i++;
                continue;
            }
            tofree = strdup(line);
            k = 0;
            while((token = strsep(&tofree, s)) != nullptr){
                if (k < trashCols){
                    if (k == 2){
                        if (i==1){
                            k++;
                            continue;
                        } else {
                            // Fetch time column
                            time = strsep(&token, " ");
                            time = token;
                            hours = strsep(&time, ":");
                            minutes = time;
                            h = strtod(hours, nullptr);
                            m = strtod(minutes, nullptr);
                            gsl_vector_set(t, i-trashRows-1, h*3600+m*60);
                            k++;
                            continue;
                        }
                    } else {
                        k++;
                        continue;
                    }
                }
                switch (i) {
                    case 1: {
                        /* Given sizes are diameter; we need radii. Also they should be microns */
                        gsl_vector_set(r, k - trashCols, strtod(token, nullptr) / 2 * 1.e-6);
                        break;
                    }
                    default: {
                        gsl_matrix_set(fv, i - trashRows - 1, k - trashCols, strtod(token, nullptr));
                        break;
                    }
                }
                k++;
            }
            i++;
        }
        fclose(f);
        gsl_vector_add_constant(t, -gsl_vector_get(t, 0));
        for (i = 1; i < t->size; i++){
            if (gsl_vector_get(t, i) - gsl_vector_get(t, i-1) == 0){
                gsl_vector_set(t, i, gsl_vector_get(t, i) + 30);
                // t_data[i] += 30;
            }
        }
        gsl_vector_set(t, 0, 0);
        if (gsl_vector_get(t, 1) - gsl_vector_get(t, 0) == 0){
            gsl_vector_set(t, 1, gsl_vector_get(t, 1) + 30);
        }
    } else perror(filename);
}

void rescaleInitial(gsl_matrix *fv, gsl_vector *r, realtype phi){
    /* Takes distributions, corresponding radii and phase fraction and rescales the distributions
     * Input args:
     *  fv  :: matrix (size MxN) of measured distributions. Distributions for given time stored in rows
     *  r   :: vector (size N) of measured droplet size classes
     *  phi :: measured phase fraction of oil in water
     *
     * fv = phi/I * f0
     * Approximate I by trapezoids: I = sum([r(i+1)-r(i)] * [f(i+1)+f(i)]) from i=0 to N-1
     */
    gsl_vector_view rowJ;
    size_t i, j;
    realtype I, ri, rii, fji, fjii;
    size_t M, N;
    M = fv->size1;
    N = fv->size2;

    for (j = 0; j < M; j++){ /* Loop over rows */
        I = 0;
        for (i = 0; i < N-1; i++){ /* Loop over columns */
            ri = gsl_vector_get(r, i); rii = gsl_vector_get(r, i+1);
            fji = gsl_matrix_get(fv, j, i); fjii = gsl_matrix_get(fv, j, i+1);
            I += ( rii-ri ) * ( fjii+fji );
        }
        I = I/2;
        rowJ = gsl_matrix_row(fv, j);
        gsl_vector_scale(&rowJ.vector, phi/I);
    }
}