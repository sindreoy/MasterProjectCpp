//
// Created by Sindre Bakke Øyen on 06.03.2018.
//

#include "Grid.h"
#include <iostream>

/* Constructors */
Grid::Grid() = default;

Grid::Grid(const Grid &g):N(g.getN()), x0(g.getX0()), x1(g.getX1()), alpha(g.getAlpha()),
                          beta(g.getBeta()), mu0(g.getMu0()),
                          xi(gsl_vector_alloc(g.getN())),
                          w(gsl_vector_alloc(g.getN())),
                          D(gsl_matrix_alloc(g.getN(), g.getN())),
                          xipBB(gsl_matrix_alloc(g.getN(), g.getN())),
                          xipBC(gsl_matrix_alloc(g.getN(), g.getN())),
                          xippBC(gsl_matrix_alloc(g.getN(), g.getN()))
{
    /* Memory has been allocated, copy values (not pointers) */
    gsl_vector_memcpy(this->xi, g.getXi());
    gsl_vector_memcpy(this->w, g.getW());
    gsl_matrix_memcpy(this->D, g.getD());
    gsl_matrix_memcpy(this->xipBB, g.getXipBB());
    gsl_matrix_memcpy(this->xipBC, g.getXipBC());
    gsl_matrix_memcpy(this->xippBC, g.getXippBC());

}

Grid::Grid(size_t N, realtype x0, realtype x1, realtype alpha, realtype beta, realtype mu0)
        : N(N), x0(x0), x1(x1), alpha(alpha), beta(beta), mu0(mu0) {
    this->w = gsl_vector_alloc(N);
    this->xi = gsl_vector_alloc(N);
    this->xipBB = gsl_matrix_alloc(N, N);
    this->xipBC = gsl_matrix_alloc(N, N);
    this->xippBC = gsl_matrix_alloc(N, N);
    this->D = gsl_matrix_alloc(N, N);

    /* Set xi, w, remap to [x0, x1] and derivative */
    this->setQuadratureRule();  /* Sets xi and w */
    this->remapGrid();                      /* Remaps from [-1, 1] to [x0, x1] */
    this->setLagrangeDerivativeMatrix();    /* In domain [x0, x1] as set by constructor */

    /* Set rescaled xis */
    this->setInterpolatedXis();
}

void Grid::coefs(size_t j, realtype *r){
    /* Returns coefficients for three-term recurrence relationship for Jacobi polynomials
 * The coefficients are
 *      r[0] = aj
 *      r[1] = bj
 *      r[2] = cj
 */
    /* Note: j runs from 0 */
    realtype a, div2, div3;
    if ((alpha == beta) == -0.5){
        r[0] = 2;
        r[1] = 0;
        r[2] = 1;
    } else {
        a = (2*j+alpha+beta);
        r[0] = (a+1)*(a+2)/(2*(j+1)*(j+alpha+beta+1));
        div2 = (2*(j+1)*(j+alpha+beta+1)*a);
        if (div2 == 0){
            r[1] = 0;
        } else{
            r[1] = (a+1)*(SUNRpowerI(alpha,2)-SUNRpowerI(beta,2))/div2;
        }
        div3 = ((j+1)*(j+alpha+beta+1)*a);
        if (div3 == 0){
            r[2] = 0;
        } else{
            r[2] = (j+alpha)*(j+beta)*(a+2)/div3;
        }
    }

}

/* Setter methods */
void Grid::setQuadratureRule() {
    /* Computes the N-point Gauss Lobatto quadrature rule with Jacobi polynomials
     * The function uses the Golub Welsch algorithm,
     *      xi = eigenvalues of Jtilde, wi = mu0 * (eigenvector vi of Jtilde)^2
     * Input args:
     * alpha:: Coefficient to determine spacing of quadrature points (alpha=beta=0 means Legendre)
     * beta :: Coefficient to determine spacing of quadrature points
     * mu0  :: Integral of weight function from a to b (a=-1, b=1 for Legendre)
     * */
    /* Declare variables */

    gsl_matrix *J, *Jtilde; /* Matrices used to obtain xi and w                 */
    gsl_matrix *I;          /* Identity matrix (size NxN)                       */
    gsl_matrix *JmI, *JpI;  /* Matrices J-I and J+I                             */
    gsl_matrix *T;          /* Matrix to solve for gammap and taup+1            */
    gsl_matrix_view subJt;  /* Submatrix of Jtilde                              */
    gsl_matrix *eVecs;      /* Matrix of eigenvectors of Jtilde                 */

    gsl_permutation *P1;    /* Permutation matrix for LU factorization          */
    gsl_permutation *P2;    /* Permutation matrix for LU factorization          */

    gsl_vector *ep;         /* Basis vector (zero except for pth element        */
    gsl_vector *eta, *mu;   /* Eigenvectors of J+-I, lambdas = +-1              */
    gsl_vector *rhs;        /* Right hand side for retrieving gammap and taup+1 */
    gsl_vector *taugamma;   /* Vector to hold gammaP and tauN+1                 */
    gsl_vector *eVals;      /* Vector of eigenvalues of Jtilde                  */

    size_t i, j;            /* Iterators                                        */
    int k = 2;              /* Signum of permutation                            */
    realtype gammai;        /* Used to store value for super and subdiagonal    */
    realtype taui;          /* Used to store value for main diagonal            */
    realtype gammaP;        /* Last value of gamma                              */
    realtype tauPp1;        /* Last value of tau                                */
    realtype etaP, muP;     /* Last values of vectors eta and mu                */
    realtype r1[3] = {0};   /* Holds coefficients from coefs function           */
    realtype r2[3] = {0};   /* Holds coefficients from coefs function           */
    realtype tmp;           /* Temporary variable                               */

    /* Allocate memory for all needed matrices and vectors */
    J       = gsl_matrix_calloc(N-1, N-1);
    Jtilde  = gsl_matrix_calloc(N, N);
    I       = gsl_matrix_alloc(N-1, N-1);
    JmI     = gsl_matrix_alloc(N-1, N-1);
    JpI     = gsl_matrix_alloc(N-1, N-1);
    T       = gsl_matrix_alloc(2, 2);
    subJt   = gsl_matrix_submatrix(Jtilde, 0, 0, J->size1, J->size2);
    eVecs   = gsl_matrix_alloc(N, N);

    P1 = gsl_permutation_alloc(N-1);
    P2 = gsl_permutation_alloc(2);

    ep      = gsl_vector_alloc(N-1);
    eta     = gsl_vector_alloc(N-1);
    mu      = gsl_vector_alloc(N-1);
    rhs     = gsl_vector_alloc(2);
    taugamma= gsl_vector_alloc(2);
    eVals   = gsl_vector_alloc(N);

    /* Construct J (normally used to calculate the Gauss quadrature rule) */
    for (i = 0; i < N-2; i++){
        this->coefs(i, r1);
        this->coefs(i+1, r2);
        taui   = r1[1]/r1[0];
        gammai = sqrt(r2[2]/(r1[0]*r2[0]));
        gsl_matrix_set(J, i, i, taui);      // Diagonal
        gsl_matrix_set(J, i, i+1, gammai);  // Superdiagonal
        gsl_matrix_set(J, i+1, i, gammai);  // Subdiagonal
    }

    /* Create all vectors and matrices to create Jtilde */
    gsl_vector_set_basis(ep, ep->size-1); // (MxN)(Nx1) = (Mx1), (J-I)eta=ep => ep in R^(Mx1)
    gsl_matrix_set_identity(I);
    gsl_matrix_memcpy(JmI, J);
    gsl_matrix_memcpy(JpI, J);
    gsl_matrix_sub(JmI, I);
    gsl_matrix_add(JpI, I);

    gsl_vector_set(rhs, 0, -1);
    gsl_vector_set(rhs, 1, 1);

    gsl_linalg_LU_decomp(JmI, P1, &k);
    gsl_linalg_LU_solve(JmI, P1, ep, mu);
    gsl_linalg_LU_decomp(JpI, P1, &k);
    gsl_linalg_LU_solve(JpI, P1, ep, eta);

    etaP = gsl_vector_get(eta, eta->size-1);
    muP = gsl_vector_get(mu, mu->size-1);

    gsl_matrix_set(T, 0, 0, 1);
    gsl_matrix_set(T, 0, 1, -etaP);
    gsl_matrix_set(T, 1, 0, 1);
    gsl_matrix_set(T, 1, 1, -muP);

    /* Enforce xi0 = -1 and xiN = 1 */
    gsl_linalg_LU_decomp(T, P2, &k);
    gsl_linalg_LU_solve(T, P2, rhs, taugamma);

    tauPp1 = gsl_vector_get(taugamma, 0);
    gammaP = sqrt(gsl_vector_get(taugamma, 1));

    gsl_matrix_swap(&subJt.matrix, J);

    gsl_matrix_set(Jtilde, Jtilde->size1-2, Jtilde->size2-1, gammaP);
    gsl_matrix_set(Jtilde, Jtilde->size1-1, Jtilde->size2-2, gammaP);
    gsl_matrix_set(Jtilde, Jtilde->size1-1, Jtilde->size2-1, tauPp1);

    /* The eigenvalues of Jtilde are the quadrature points and
     * the first value of each eigenvector is used to obtain the quadrature weights
     * */
    gsl_eigen_symmv_workspace *ws;
    ws = gsl_eigen_symmv_alloc(N);
    gsl_eigen_symmv(Jtilde, eVals, eVecs, ws);

    /* Find quadrature points */
    gsl_vector_swap(xi, eVals);
    /* Find weights */
    for (i=0; i < N; i++){
        tmp = gsl_matrix_get(eVecs, 0, i);
        tmp *= tmp;
        gsl_vector_set(w, i, tmp*mu0);
    }

    /* xi and w are reversed. 3/4 of the xi and w are sorted. Sort with insertion sort */
    gsl_vector_reverse(xi);
    gsl_vector_reverse(w);
    i = 1;
    while (i < xi->size){
        j = i;
        while ((j > 0) && (gsl_vector_get(xi, j-1) > gsl_vector_get(xi, j))){
            /* Some value is less than previous: swap elements until it is not */
            gsl_vector_swap_elements(xi, j, j-1);
            gsl_vector_swap_elements(w, j, j-1);  /* Remember that w_i corresponds to xi_i: swap accordingly */
            j--;
        }
        i++;
    }

    /* END OF FUNCTION: FREE ALLOCATED MEMORY */
    gsl_matrix_free(J);
    gsl_matrix_free(Jtilde);
    gsl_matrix_free(I);
    gsl_matrix_free(JmI);
    gsl_matrix_free(JpI);
    gsl_matrix_free(T);
    gsl_matrix_free(eVecs);

    gsl_permutation_free(P1);
    gsl_permutation_free(P2);

    gsl_vector_free(ep);
    gsl_vector_free(eta);
    gsl_vector_free(mu);
    gsl_vector_free(rhs);
    gsl_vector_free(taugamma);
    gsl_vector_free(eVals);
    gsl_eigen_symmv_free(ws);
}

void Grid::remapGrid(){
    /* Remap xi from [-1,1] to physical domain [this->x0, this->x1] */
    realtype a0 = gsl_vector_get(this->xi, 0);
    realtype b0 = gsl_vector_get(this->xi, this->N-1);

    gsl_vector_scale(this->xi, (this->x1-this->x0)/(b0-a0));
    gsl_vector_add_constant(this->xi, -a0*(this->x1-this->x0)/(b0-a0));

    gsl_vector_scale(this->w, (this->x1-this->x0)/(b0-a0));
}

void Grid::setLagrangeDerivativeMatrix(){
    realtype s = 1;
    realtype wi, wj;
    size_t i, j, k;
    for (i = 1; i < this->N+1; i++){       /* Loop on rows */
        s = 1;
        for (k = 0; k < this->N; k++){
            if (k != i-1){
                s /= gsl_vector_get(this->xi, i-1) - gsl_vector_get(this->xi, k);
            }
        }
        wi = s;
        for (j = 1; j < this->N+1; j++){   /* Loop on columns (l_j(x_i)) */
            if (i==j){ /* Diagonal */
                s = 0;
                for (k = 0; k < this->N; k++){
                    if (k != i-1){
                        s += 1 / (gsl_vector_get(this->xi, i-1) - gsl_vector_get(this->xi, k));
                    }
                }
                gsl_matrix_set(D, i-1, j-1, s);
            } else{ /* Off-diagonal */
                s = 1;
                for (k = 0; k < this->N; k++){
                    if (k != j-1){
                        s /= gsl_vector_get(this->xi, j-1) - gsl_vector_get(this->xi, k);
                    }
                }
                wj = s;
                gsl_matrix_set(D, j-1, i-1, wi /
                        (wj * (gsl_vector_get(this->xi, j-1)-gsl_vector_get(this->xi, i-1))));
            }

        }
    }
}

void Grid::setInterpolatedXis(){
    size_t i, j;
    realtype ai, bj, tmp;
    realtype exp1, exp2;
    for (i = 0; i < this->N; i++){
        ai = gsl_vector_get(xi, i);
        for (j = 0; j < this->N; j++){
            bj   = gsl_vector_get(xi, j);
            exp1 = SUNRpowerR(bj, 3.0)/2;
            tmp  = 1.0 - exp1;
            exp2 = SUNRpowerR(tmp, 1.0/3);
            gsl_matrix_set(xipBB, i, j, (1.0-ai)*bj+ai);                    /* Breakage birth    */
            gsl_matrix_set(xipBC, i, j, SUNRpowerR(2.0, -1.0/3)*ai*bj);     /* Coalescence birth */
            gsl_matrix_set(xippBC, i, j, ai*exp2);                          /* Coalescence birth */
        }
    }
}

/* Getter methods */
size_t Grid::getN() const {
    return N;
}

realtype Grid::getX0() const {
    return x0;
}

realtype Grid::getX1() const {
    return x1;
}

realtype Grid::getAlpha() const {
    return alpha;
}

realtype Grid::getBeta() const {
    return beta;
}

realtype Grid::getMu0() const {
    return mu0;
}

gsl_vector *Grid::getXi() const {
    return xi;
}

gsl_vector *Grid::getW() const {
    return w;
}

gsl_matrix *Grid::getD() const {
    return D;
}

gsl_matrix *Grid::getXipBB() const {
    return xipBB;
}

gsl_matrix *Grid::getXipBC() const {
    return xipBC;
}

gsl_matrix *Grid::getXippBC() const {
    return xippBC;
}

std::ostream &operator<<(std::ostream &os, const Grid &grid) {
    size_t i, j, N = grid.getN();
    gsl_vector *xi = grid.getXi();
    gsl_vector *w  = grid.getW();
    gsl_matrix *D  = grid.getD();
    gsl_matrix *xipBB = grid.getXipBB();
    gsl_matrix *xipBC = grid.getXipBC();
    gsl_matrix *xippBC = grid.getXippBC();

    os << "xi (Quadrature points):\n";
    for (i = 0; i < N; i++){
        os << std::setprecision(3) << gsl_vector_get(xi, i) << "\t";
    }
    os << "\n\nw (Quadrature weights):\n";
    for (i = 0; i < N; i++){
        os << std::setprecision(3) << gsl_vector_get(w, i) << "\t";
    }
    os << "\n\nD (Lagrange derivative matrix):\n";
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            os << std::setw(8) << std::setprecision(3) << gsl_matrix_get(D, i, j) << "\t";
        }
        os << "\n";
    }
    os << "\n\nxipBB (Quadrature points birth breakage):\n";
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            os << std::setw(8) << std::setprecision(3) << gsl_matrix_get(xipBB, i, j) << "\t";
        }
        os << "\n";
    }
    os << "\n\nxipBC (Quadrature points birth coalescence):\n";
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            os << std::setw(8) << std::setprecision(3) << gsl_matrix_get(xipBC, i, j) << "\t";
        }
        os << "\n";
    }
    os << "\n\nxippBC (Quadrature points birth coalescence):\n";
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            os << std::setw(8) << std::setprecision(3) << gsl_matrix_get(xippBC, i, j) << "\t";
        }
        os << "\n";
    }
    return os;
}

/* Destructors */
Grid::~Grid(){
    gsl_vector_free(this->xi);
    gsl_vector_free(this->w);
    gsl_matrix_free(this->D);
    gsl_matrix_free(this->xipBB);
    gsl_matrix_free(this->xipBC);
    gsl_matrix_free(this->xippBC);
}
