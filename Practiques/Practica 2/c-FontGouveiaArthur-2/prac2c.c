/* Font Gouveia Arthur 30306373D */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prac2funs.h"

#define PI 3.14159265358979323846

int main(void) {
    int iterMax, res, i, countSol = 0;
    double *sol, tol;
    double **sols, *Q , a, b;
    
    /* Memory Allocation */
    sols = (double**) malloc(3*sizeof(double));
    sol = (double*) malloc(3*sizeof(double));
    Q = (double*) malloc(3*sizeof(double));

    for (i = 0; i < 3; i++){
        sols[i] = (double *) malloc(3*sizeof(double) );
    }
    
    /* Por ejemplo, tol = 1e-4 */
    printf("Doneu la tolerancia: ");
    scanf("%le", &tol);
 
    /* Por ejemplo, iterMax = 10 */
    printf("Doneu el nombre màxim d'iteracions: ");
    scanf("%d", &iterMax);

    printf("Doneu el rang [a,b]:\n");
    scanf("%le", &a);
    scanf("%le", &b);
    
    srand48(0);	
    /* Buscamos possibles raizes de la funcion f(x) */
    for (i = 0; i < 100; i++) {
	    Q[0] = (drand48() * fabs(b-a)) + a;
            Q[1] = (drand48() * fabs(b-a)) + a;
	    Q[2] = (drand48() * fabs(b-a)) + a;
	    /* Aplicamos el método de Newton en ese sub intervalo */
	    res = newton3 (Q , sol, tol, iterMax);
	    if (res == 0) {
		/* Guardamos la solucion el en vector de solucions*/
	        sols[countSol] = sol;
                countSol++;
	    }
	    else {
	    	printf("No s’ha trobat una aproximació prou bona del zero.");
    	    }
    }
	
    return 0;
}


int newton3(double *x, double *sol, double tol, int iter){
    int i, j, k, a;
    double *x0,*x1;
    double **Df, **DfT, **C;
    double *deltaXn, *f, *b;
    
    /* Memory Allocation */
    x0 = (double*) malloc(3*sizeof(double));
    x1 = (double*) malloc(3*sizeof(double));
    Df = (double **) malloc(3*sizeof(double *) );
    DfT = (double **) malloc(3*sizeof(double *) );
    C = (double **) malloc(3*sizeof(double *) );
    deltaXn = (double *) malloc(3*sizeof(double *) );
    f = (double *) malloc(3*sizeof(double *) );
    b = (double *) malloc(3*sizeof(double *) );

    if (x0 == NULL || x1 == NULL || Df == NULL || DfT == NULL || 
         b == NULL || C == NULL || deltaXn == NULL || f == NULL) {
        printf("No hi ha prou memoria");
        exit (1);
    }

    for (i = 0; i < 3; i++){
        Df[i] = (double *) malloc(3*sizeof(double) );
        DfT[i] = (double *) malloc(3*sizeof(double) );
        C[i] = (double *) malloc(3*sizeof(double *) );
        if ( Df[i] == NULL || DfT[i] == NULL){
            printf("No hi ha prou memoria");
            exit (2);
        }
    }

    /* Metodo Newton3 */    

    /* Asignamos la aproximación inicial x0 = x */
    for (i = 0; i < 3; i++){
   	x0[i] = x[i];
    }

    for(k = 0; k < iter; k++){
	/* Calculamos deltaXn resolvendo un sistema lineal */
	
	F(x,f);	    
	dF(x,Df);
        /* Transposem la matriu Df */
	for (i = 0; i < 3; i++){
	    for (j = 0;j < 3;j++){
	        DfT[j][i] = Df[i][j];
	    }
	}
        prodMatMat(3,3,3,DfT,Df,C);
        prodMatVec(3,3,DfT,f,b);
        /* Apliquem la Descompoició LDL^t */  
        a = ldlt(3,C,0.1);
	    if (a == 1) {
		printf("No es diagonalizable");
	    }else {
		for (i = 0; i < 3; i++){
		    for (j = 0; j < 3; j++){
		        printf("%16.7e ", C[i][j]);
		    }
		    printf("\n");
		}
		
		/* Resoldrem el sisteme triangular inferior */  
		resTinf(3, C, b, deltaXn);
		printf("\n Vector solució deltaXn = \n");
		for (i = 0; i < 3; i++){
		    printf( "%16.7e ", deltaXn[i]);
		}
		printf("\n");
	}

	/* Aqui ya tenemos deltaXn */
      	/* Xn+1 = Xn - deltaXn */
        for (i = 0; i < 3; i++){
   	    x1[i] = x0[i] - deltaXn[i];
        }
      	
      	/* Si la norma2 (3, deltaXn) < tol, se ha encontrado una 			aproximación del cero suficientemente buena */
      	if (norma2(3, deltaXn) < tol ) {
      	    /* Asignamos sol = Xn+1 y la función retorna 0*/
      	    *sol = *x1;
      	    return 0;
      	}
      	/* Se atualiza la Xn */
	for (i = 0; i < 3; i++){
   	    x0[i] = x1[i];
        }
    }
    /* No se ha encontrado una aproximación del cero suficientemente 	 	buena */
    return 1;
}

/* Sistema : x + y + z = 1
	         y + z = 0 
            x² + 0.75y = 0 */
void F(double *x, double *f) {
    f[0] = x[0] + x[1] + x[2] - 1;
    f[1] = x[1] + x[2];
    f[2] = pow(x[0], 2) + x[1] * 0.75;
}

/* Matriz Df de derivadas parciales */
void dF(double *x, double **df) {
   df[0][0] = 1; 
   df[0][1] = 1;
   df[0][2] = 1;
   df[1][0] = 0;
   df[0][1] = 1;
   df[0][2] = 1;
   df[2][0] = 2* x[0];
   df[2][1] = 0.75;
   df[2][1] = 0;
}
