/* Font Gouveia Arthur 30306373D */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prac1funs.h"

int main(void){
    int n, i, j;
    double **L, **T, *x, *b, *vr, total = 0;
    printf("Doneu les dimensions de la matriu (n,n) = \n");
    scanf("%d", &n);
    L = (double **) malloc(n*sizeof(double *) );
    T = (double **) malloc(n*sizeof(double *) );
    if (L == NULL || T == NULL){
        printf("No hi ha prou memoria");
        exit (1);
    }
    for (i = 0; i < n; i++){
        L[i] = (double *) malloc(n*sizeof(double) );
        T[i] = (double *) malloc(n*sizeof(double) );
        if ( L[i] == NULL || T[i] == NULL){
            printf("No hi ha prou memoria");
            exit (2);
        }
    }
    x = (double *) malloc(n*sizeof(double) );
    b = (double *) malloc(n*sizeof(double) );
    vr = (double *) malloc(n*sizeof(double) );
    if (x == NULL || b == NULL || vr == NULL){
        printf("No hi ha prou memoria");
        exit (3);
    }
    printf("Doneu els (%d x %d) elements de la matriu L \n", n, n);
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            scanf("%le", &L[i][j]);
        }
        if (L[i][i] != 1){
            printf("No te diagonal de 1's");
            exit (4);
        }
    }
    
    printf("Doneu els %d elements del vector b \n", n);
    for (i = 0; i < n; i++)
        scanf("%le", &b[i]);

    /* Resoldrem el sistema triangular inferior */ 
    resTinf(n, L, b, x);
    
    /* Vector solució */   
    printf("\n Vector solució x = \n");
    for (i=0;i<n;i++)
        printf( "%16.7e ", x[i]);


    /* Valor residual */   
    printf("\n Valor residual ||Lx - b||^2 = \n");

    /* Producte Matriu per Vector */ 
    prodMatVec (n,n,L,x,vr);
    for (i = 0; i < n; i++){
        vr[i] = vr[i] - b[i];
        vr[i] = pow(vr[i],2);
        total += vr[i];
    }
    
    total = sqrt(total);
    printf( "%16.7e ", total);
    
    /* Transposem la matriu L */
    for (i=0; i<n; i++){
        for (j=0;j<n;j++){
            T[j][i] = L[i][j];
        }
    }
    
    printf("\n Matriu L transposada = \n");
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            printf(" %16.7e ", T[i][j]);
        }
        printf("\n");
    }
    
    /* Resoldrem el sistema triangular superior */ 
    resTsup(n,T,b,x);
    
    /* Vector solució de la transposada*/  
    printf("\n Vector solució x de la Transposada = \n");
    for (i = 0; i < n; i++){
        printf( "%16.7e ", x[i]);
    }
    printf("\n");

    for (i = 0; i < n; i++) {
        free (L[i]);
        free (T[i]);
    }   
    free(L);
    free(T);
    free(x);
    free(b);
    free(vr);
    return 0;
}
