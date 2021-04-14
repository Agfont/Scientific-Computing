/* Font Gouveia Arthur 30306373D */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prac2funs.h"

#define PI 3.14159265358979323846

int main(void) {
    int iterMax, res, i, countSol = 0;
    double sol, aproxInicial, tol;
    double *sols, n, M, h, a, b, c, prev, current;
    
    /* Por ejemplo, tol = 1e-4 */
    printf("Doneu la tolerancia: ");
    scanf("%le", &tol);

    /* Por ejemplo, iterMax = 10 */
    printf("Doneu el nombre màxim d'iteracions: ");
    scanf("%d", &iterMax);

    /* El nº máximo de raízes reales és el grado del polinomio */    
    printf("Doneu el grau del polinomi: ");
    scanf("%le", &n);
    
    /* Depende de cada caso */
    printf("Doneu la constant M: ");
    scanf("%le", &M);
  
    /* Memory Allocation */
    sols = (double*) malloc(n*sizeof(double));

    /* Definimos el rango, el paso y siguiente */
    a = -M;
    b = M;
    h = 2*M/1000;
    c = a + h;
    prev = f(a);
    
    /* Si a és raíz de la funcion f(x) */
    if (prev == 0) {
        sols[countSol] = prev;
        countSol++;
    }
    
    /* Buscamos possibles raizes de la funcion f(x) */
    while (c < b && countSol < n) {
	current = f(c);
	/* Miramos si existe una solución en el sub intervalo */
	if (prev * current < 0) {
	    aproxInicial = 0.5 * (a + c);
	    /* Aplicamos el método de Newton en ese sub intervalo */
	    res = newton (aproxInicial , &sol, tol, iterMax);
	    if (res == 0) {
		/* Guardamos la solucion el en vector de solucions*/
	        sols[countSol] = sol;
                countSol++;
	    }
	    else {
	    	printf("No s’ha trobat una aproximació prou bona del zero.");
    	    }
	}
	/* Atualizamos las variables */
	prev = current;
	a = c;
	c = c + h;
    }
    
    /* Mostramos las soluciones */
    printf("Solucions: \n");
    for (i = 0; i < countSol; i++){
        printf("%16.7e\n", sols[i]);
    }
    
    free(sols);

    return 0;
}


int newton(double x, double *sol, double tol, int iter){
    int i;
    double x0,x1;
    /* Asignamos la aproximación inicial x0 = x */
    x0 = x;

    for(i = 0; i < iter; i++){
      	/* Xn+1 = Xn - F(Xn) / dF(Xn) */
        x1 = x0 - f(x0) / df(x0);
      	
      	/* Si |Xn+1 - Xn| < tol o |F(Xn)| < tol, se ha encontrado una aproximación del cero suficientemente buena */
      	if (fabs(x1-x0) < tol || fabs(f(x0)) < tol) {
      	    /* Asignamos sol = Xn+1 y la función retorna 0*/
      	    *sol = x1;
      	    return 0;
      	}
      	/* Si |dF(Xn)| < tol, no se ha encontrado una aproximación del cero suficientemente buena */
      	if (fabs(df(x0)) < tol) {
      	    /* La función retorna 1 */
      	    return 1;
      	}
      
      	/* Se atualiza la Xn */
      	x0 = x1;
    }
    /* No se ha encontrado una aproximación del cero suficientemente 	 	buena */
    return 1;
}

/* P(x) = x² - 1 ; M = 2*/
double f(double x) {
    return pow(x, 2) - 1;
}

/* P'(x) = 2x */
double df(double x) {
    return 2*x;
}

/* Otros polinomios */

/* b) -------- P(x) = x^3 − x ; M = 2 ---------------- */
/* double f(double x) {
    return pow(x,3) - x
}

double df(double x) {
    return 3 * pow(x, 2) - 1;
}*/

/* c) -------- P(x) = 3x^3 - x + 1 ; M = 5/3 -------------*/
/* double f(double x) {
    return 3 * pow(x,3) - x + 1
}

double df(double x) {
    return 6 * pow(x,2) - 1
}*/

/* d) -------- P(x) = x^4 + 1 ; M = 2 ------------*/
/* double f(double x) {
    return pow(x,4) + 1
}

double df(double x) {
    return 4 * pow(x,3)
}*/

/* e) -------- P(x) = Prod{i=1, i<=6}(x - 10i/i+1) ; M = 1250870 -- */
/* f) -------- P(x) = Prod{i=1, i<=6}(x - i+1/10*i) ; M = 8.756 --- */
