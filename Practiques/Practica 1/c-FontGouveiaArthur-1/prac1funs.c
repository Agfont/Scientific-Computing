/* Font Gouveia Arthur 30306373D */

#include <math.h>
#include "prac1funs.h"

void resTinf (int n,double **L, double *b, double *x) {
    int i, k;

    /* Sucessión X[i] = B[i] - sum(L[i][k] * x[k]) */
    for (i = 0; i < n; i++) {
        x[i] = b[i];
        
        for (k = 0; k < i; k++) {
            x[i] -= (L[i][k] * x[k]);
        }
    }
}

void resTsup (int n,double **U, double *b, double *x) {
    int i, k;

    /* Sucessión X[i] = B[i] - sum(U[i][k] * x[k]) */
    for (i = n-1; i >= 0; i--) {
        x[i] = b[i];
            
	for (k = n-1; k > i; k--) {
	    x[i] -= (U[i][k] * x[k]);
	}
    }
    
}
void prodMatVec (int m, int n, double **A, double *x, double *y) {
    int i, j;
    /* Ponemos el vector solución a 0 */
    for (i = 0; i < n; i++) {
 	y[i] = 0.f;
        /* Hacemos el producto de matrices "linea*coluna" */
	for (j = 0; j < n; j++) {
	    y[i] += A[i][j] * x[j];
	}
    }
}
void prodMatMat (int m, int n, int p, double **A, double **B, double **C) {
    int i, j, k;
    /* Ponemos la matriz solución a 0 */
    for (i = 0; i < m; i++) {
	for (j = 0; j < p; j++) {
  	    C[i][j] = 0.f;
	    for(k = 0; k < n; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
double norma2 (int n, double *z) {

    /* Declaración de variables */
   int i; 
   double norma = 0;
    
    /* Norma = Raíz de la suma de los cuadrados */
    for (i = 0; i < n; i++) {
	norma += pow(z[i], 2);
    }
    norma = sqrt(norma);
    return norma;
}
