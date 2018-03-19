//
// Created by Sindre Bakke Øyen on 18.03.2018.
//

#include <gsl/gsl_matrix.h>
#include "PBModel.h"

PBModel::PBModel() = default;

//PBModel::PBModel(const gsl_matrix &fv0, const gsl_vector &psi_vec, const gsl_vector &r, const gsl_vector &t,
//             const Grid &g, const Kernels &k, const Constants &c):
//psi_mat(SUNDenseMatrix(fv0.size1, fv0.size2)), psi(N_VNew_Serial(psi_vec.size)), r(N_VNew_Serial(r.size)),
//time(N_VNew_Serial(t.size)), grid(g), kerns(k), consts(c), LS(nullptr), cvode_mem(nullptr)
//{
//    size_t i, j, N;
//    N = g.getN();
//
//    /* Copy elements from fv0 [MxN] into psi_mat */
//    for (i = 0; i < N; i++){
//
//    }
//}
PBModel::PBModel(char const *f, const Grid &g, const Kernels &k,
             const Constants &c, const SystemProperties &s):
        filename(f), grid(g), kerns(k), consts(c), sysProps(s){
    /* Get rows and columns of csv file and initialize M and N respectively */
    size_t rows=0, cols=0;
    this->getRowsAndCols(rows, cols);
    this->M      = rows;
    this->N      = cols;

    /* Allocate memory for member variables */
    this->t      = gsl_vector_alloc(this->M);
    this->r      = gsl_vector_alloc(this->N);
    this->fv     = gsl_matrix_alloc(this->M, this->N);
    this->gslPsi = gsl_matrix_alloc(this->M, this->N);

    /* Set t, r, fv and rescale fv to phase fraction */
    this->getDistributions(rows, cols);
    this->rescaleInitial();

    /* Nondimensionalize fv (psi = fv*Rmax) */
//    gsl_vector_view fv0 = gsl_matrix_row(this->fv, 0);
//    gsl_matrix_scale(this->gslPsi, sysProps.getRm());

}

void PBModel::getRowsAndCols(size_t &rows, size_t &cols){
    FILE *f = fopen(filename, "r");
    if (f != nullptr) {
        // Initialize variables
        rows = 0;
        cols = 0;
        size_t i=0, j=1000, tmp = 0;
        bool flag = false;
        realtype val = 0;
        const char s[2] = ",";
        const int bufSize = 1000;
        char line[bufSize], *toFree;

        // Now loop over lines
        while(fgets(line, sizeof line, f) != nullptr) {
            rows++;
            tmp = 0;
            i = 0;
            toFree = strdup(line);
            while ((strsep(&toFree, s)) != nullptr) {
                tmp++;
                if (rows > TRASHROWS && tmp > TRASHCOLS){
                    if (toFree != nullptr) {
                        val = strtod(toFree, nullptr);
                    } else val = 1;
                    if (val != 0) flag = true;
                    else if (flag && val == 0){
                        i++;
                    }
                }
            }
            flag = false;
            if (rows > TRASHROWS && j > i) j = i;
            if (tmp > cols) {
                cols = tmp;
            }
        }
        if (j > TRUNCATETHRESHOLD){
            cols = cols - TRASHCOLS - j + TRUNCATETHRESHOLD;
        }
        rows = rows - TRASHROWS;
        fclose(f);
    } else {
        perror(filename);
    }
}

void PBModel::getDistributions(size_t &rows, size_t &cols) {
    /* Declare needed variables */
    char *hours, *minutes;
    realtype h, m;
    size_t i = 0, k = 0;; // Index variables

    const char s[2] = ",";
    const size_t bufSize = 10*(cols);
    char line[bufSize], *token, *toFree, *timeStr;

    /* Open file and start reading */
    FILE *f = fopen(filename, "r");
    if (f != nullptr) {
        while (fgets(line, sizeof line, f) != nullptr) {
            if (i == 0) {
                i++;
                continue;
            }
            toFree = strdup(line);
            k = 0;
            while((token = strsep(&toFree, s)) != nullptr){
                if (k < TRASHCOLS){
                    if (k == 2){
                        if (i==1){
                            k++;
                            continue;
                        } else {
                            // Fetch time column
                            timeStr = strsep(&token, " ");
                            timeStr = token;
                            hours = strsep(&timeStr, ":");
                            minutes = timeStr;
                            h = strtod(hours, nullptr);
                            m = strtod(minutes, nullptr);
                            gsl_vector_set(t, i-TRASHROWS, h*3600+m*60);
                            k++;
                            continue;
                        }
                    } else {
                        k++;
                        continue;
                    }
                }
                if (k-TRASHCOLS < this->N) {
                    switch (i) {
                        case 1: {
                            /* Given sizes are diameter; we need radii. Also they should be microns */
                            gsl_vector_set(r, k - TRASHCOLS, strtod(token, nullptr) / 2 * 1.e-6);
                            break;
                        }
                        default: {
                            gsl_matrix_set(fv, i - TRASHROWS, k - TRASHCOLS, strtod(token, nullptr));
                            break;
                        }
                    }
                    k++;
                } else continue;
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

void PBModel::rescaleInitial(){
    /* Takes distributions, corresponding radii and phase fraction and rescales the distributions
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
        gsl_vector_scale(&rowJ.vector, PHI/I);
    }
}

gsl_matrix *PBModel::getFv() const {
    return fv;
}

gsl_vector *PBModel::getR() const {
    return r;
}

gsl_vector *PBModel::getT() const {
    return t;
}

size_t PBModel::getM() const {
    return M;
}

size_t PBModel::getN() const {
    return N;
}

void PBModel::printDistribution(){
    size_t i, j;
    std::cout << "The droplet size density distribution:" << std::endl;
    for (i=0;i<M;i++){
        for(j=0;j<N;j++){
            std::cout << std::setw(4) << std::setprecision(4) << gsl_matrix_get(fv, i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

void PBModel::printDimensions(){
    std::cout << "Rows: " << this->M << ", Columns: " << this->N << std::endl;
}

void PBModel::printTime() {
    size_t i = 0;
    std::cout << "Time of measurement: " << std::endl;
    for (i = 0; i < this->M; i++){
        std::cout << std::setw(5) << std::setprecision(4) << gsl_vector_get(this->t, i);
    }
    std::cout << std::endl;
}

void PBModel::printSizeClasses() {
    size_t i = 0;
    std::cout << "Measured size classes: " << std::endl;
    for (i = 0; i < this->N; i++){
        std::cout << std::setw(10) << std::setprecision(7) << gsl_vector_get(this->r, i);
    }
    std::cout << std::endl;
}

PBModel::~PBModel(){
    gsl_matrix_free(this->fv);
}
