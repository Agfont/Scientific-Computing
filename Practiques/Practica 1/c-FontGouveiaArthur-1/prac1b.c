/* Font Gouveia Arthur 30306373D */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prac1funs.h"

int ldlt(int n, double **A, double tol);

/* Descomposició LDL^t */
int ldlt(int n, double **A, double tol){
    int i, j, k;

    for (k=0; k <= n-1; k++){
	/* Sucessió A[k][k] = A[k][k] - sum(A[k][j]² * A[j][j]) */
        for (j=0; j<= k-1; j++){
            A[k][k] -= pow(A[k][j],2) * A[j][j];
        }
 
	/* Sucessió (A[i][k] - sum(A[i][j] * A[k][j]* A[j][j]) */
        for (i=k+1; i<= n-1; i++){
            for (j=0; j<= k-1; j++){
                A[i][k] -= (A[i][j] * A[k][j] * A[j][j]);
            }
            /* Sucessió A[i][k] = 1/A[k][k] */
            A[i][k] = A[i][k] * (1/A[k][k]);
            A[k][i] = A[i][k];
        }
    }
    
    /* Verifiquem se s'ha pogut descomposar la matriu */
    for (k = 0; k<n; k++){
	/* Si no s'ha pogut, retorna 1 */
        if (fabs(A[k][k]) < tol )
            return 1;
    }
    
    /* Si s'ha pogut, retorna 0 */
    return 0;
}


int main(void){
    int n, i, j, a;
    double **A, **T, *x, *y, *z, *b, total = 0;
    printf("Doneu les dimensions de la matriu (n,n) = \n");
    scanf("%d", &n);
    A = (double **) malloc(n*sizeof(double *));
    T = (double **) malloc(n*sizeof(double *));
    
    if (A == NULL){
        printf("No hi ha prou memoria");
        return 1;
    }
    for (i = 0; i < n; i++){
        A[i] = (double *) malloc(n*sizeof(double));
        T[i] = (double *) malloc(n*sizeof(double));
        if ( A[i] == NULL){
            printf("No hi ha prou memoria");
            return 2;
        }
    }
    x = (double *) malloc(n*sizeof(double));
    y = (double *) malloc(n*sizeof(double));
    z = (double *) malloc(n*sizeof(double));
    b = (double *) malloc(n*sizeof(double));
    if ( x == NULL || b == NULL || y == NULL || z == NULL){
        printf("No hi ha prou memoria");
        return 3;
    }
    printf("Doneu el (%d x %d) element de la matriu L \n", n, n);
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            scanf("%le", &A[i][j]);
        }
    }

    /* Comprovem si no és simetrica */  
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            T[i][j] = A[j][i];
	     if (T[i][j] != A[i][j]){
                printf("No és simetrica");
                return 4;
             }
        }
    }
    
    printf("Doneu els %d elements del vector b \n", n);
    for (i = 0; i < n; i++)
        scanf("%le", &b[i]);
    
    /* Apliquem la Descompoició LDL^t */  
    a = ldlt(n,A,0.1);

    if (a == 1) {
        printf("No es diagonalizable");
        for (i = 0; i < n; i++)
            free (A[i]);
        
        free(A);
	free(x);
	free(y);
        free(z);
        free(b);
        return 5;
        
    }
    else {
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                printf ("%16.7e ", A[i][j]);
            }
            printf("\n");
        }

	/* Resoldrem el sisteme triangular inferior */          
        resTinf(n, A, z, b);
        
        for (i=0; i<n; i++) {
            y[i] = z[i]/A[i][i];
        }
        
	/* Resoldrem el sisteme triangular superior */  
        resTsup(n, A, x, y);
        
        /* Valor residual */
        printf("\n Valor residual ||Ax - b||^2 = \n");
	/* Producte Matriu * Vector */
        prodMatVec (n,n,A,x,y);
        for (i=0; i <n; i++){
            y[i] = y[i] - b[i];
            y[i] = pow(y[i],2);
            total += y[i];
        }
        
        total = sqrt(total);
        printf( "%16.7e ", total);
        
        /* Vector solució */  
        printf("\n Vector solució x = \n");
        for (i=0;i<n;i++){
            printf( "%16.7e ", x[i]);
        }
    }
    printf("\n");
    
    for (i = 0; i < n; i++) {
        free (A[i]);
    }
    free(A);
    free(z);
    free(x);
    free(b);
    free(y);

    return 0;
}


