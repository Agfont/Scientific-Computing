/* Font Gouveia Arthur 30306373D */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
    int m, n, i, j, a;
    double **A, **T,**C, *x, *z, *b, *vr, *v, total = 0;
    printf("Doneu les dimensions de la matriu (m,n) = \n");
    scanf("%d %d", &m,&n);
    A = (double **) malloc( m*sizeof(double *) );
    T = (double **) malloc( n*sizeof(double *) );
    C = (double **) malloc( n*sizeof(double *) );
    if (A == NULL || T == NULL || C == NULL){
        printf("No hi ha prou memoria");
        exit (1);
    }
    for (i = 0; i < m; i++){
        A[i] = (double *) malloc( n*sizeof(double) );
    }

    for (i = 0; i < n; i++){
        T[i] = (double *) malloc( m*sizeof(double) );
        C[i] = (double *) malloc( n*sizeof(double) );
        if (A[i] == NULL || T[i] == NULL || C[i] == NULL){
            printf("No hi ha prou memoria");
            exit (2);
        }
    }
    x = (double *) malloc( n*sizeof(double) );
    b = (double *) malloc( m*sizeof(double) );
    vr = (double *) malloc( n*sizeof(double) );
    z = (double *) malloc( n*sizeof(double) );
    v = (double *) malloc( n*sizeof(double) );

    if (x == NULL || b == NULL || vr == NULL || z == NULL || v == NULL){
        printf("No hi ha prou memoria");
        exit (3);
    }
    printf("Doneu els (%d x %d) elements de la matriu A \n", m, n);
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            scanf("%le", &A[i][j]);
        }
    }
    
    printf("Doneu els %d elements del vector b \n", m);
    for (i = 0; i < m; i++)
        scanf("%le", &b[i]);
	
    /* Transposem la matriu A */
    for (i= 0; i<m; i++){
        for (j=0;j<n;j++){
            T[j][i] = A[i][j];
        }
    }
    
    prodMatMat(n,m,n,T,A,C);
    prodMatVec(n,n,C,b,v);
	
    printf("\n Vector solució x = \n");
    for (i = 0; i < n; i++){
	for (j = 0; j < n; j++) {
            printf( "%16.7e ", C[i][j]);
	}
    }

    /* Apliquem la Descompoició LDL^t */  
    a = ldlt(n,C,0.1);
    if (a == 1) {
        printf("No es diagonalizable");
	for (i = 0; i < m; i++){
            free (A[i]);
	}
        for (i = 0; i < n; i++){
            free (T[i]);
            free (C[i]);
        }
        
        free(A);
        free(T);
        free(C);
        
        free(z);
        free(x);
        free(b);
        free(vr);
        free(v);

        exit (4);
        
    } else{
        
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                printf("%16.7e ", C[i][j]);
            }
            printf("\n");
        }
	
	/* Resoldrem el sisteme triangular inferior */  
        resTinf(n, C, v, z);
    
        for (i=0; i<n; i++) {
            vr[i] = z[i]/C[i][i];
        }
    	
	/* Resoldrem el sisteme triangular superior */  
        resTsup(n, C, x, vr);
    
        /* Valor residual */
        printf("\n Valor residual ||Ax - b||^2 = \n");
	/* Producte Matriu * Vector */
        prodMatVec (n,n,C,x,vr);
        
	for (i=0; i < n; i++){
            vr[i] = vr[i] - v[i];
            vr[i] = pow(vr[i],2);
            total += vr[i];
        }
        total = sqrt(total);
        printf( "%16.7e ", total);
        
	/* Vector solució */
        printf("\n Vector solució x = \n");
        for (i = 0; i < n; i++){
            printf( "%16.7e ", x[i]);
        }
        printf("\n");
        
	for (i = 0; i < m; i++){
            free (A[i]);
        }
        for (i = 0; i < n; i++){
            free (T[i]);
            free (C[i]);
        }
        
        free(A);
        free(T);
        free(C);
        
        free(z);
        free(x);
        free(b);
        free(vr);
        free(v);

        return 0;
    }
}
