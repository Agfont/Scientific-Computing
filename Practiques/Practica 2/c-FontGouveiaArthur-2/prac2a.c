/* Font Gouveia Arthur 30306373D */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prac2funs.h"

#define PI 3.14159265358979323846

int main(void) {
    int iterMax, resultado;
    double sol, aproxInicial, tol;

    /* Por ejemplo: tol = 1e-4 */
    printf("Doneu la tolerancia: ");
    scanf("%le", &tol);

    /* Por ejemplo: iterMax = 10 */
    printf("Doneu el nombre màxim d'iteracions: ");
    scanf("%d", &iterMax);

    /* Depende de cada caso */
    printf("Doneu l'aproximació inicial: ");
    scanf("%le", &aproxInicial);

    /* Resolvemos utilizando el metodo de Newton */
    resultado = newton (aproxInicial , &sol, tol, iterMax);

    if (resultado == 0) {
        printf("Solució = ");
	printf("%16.7e", sol);
	printf("\n");
    }
    else {
	printf("No s’ha trobat una aproximació prou bona del zero.");
    }

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

/* a) f(x) = x² + sin'(x) - PI ; aproxInicial = 2 */
double f(double x) {
    return pow(x, 2) + sin(x) - PI;
}

/* f'(x) = 2x + cos(x)*/
double df(double x) {
    return (2 * x) + cos(x);
}

/* Otras funciones */

/* b) -------- f (x) = 1 − log(x) ; aproxInicial = 3 --------*/
/* double f(double x) {
    return 1 - log(x)
}

double df(double x) {
    return -1/x;
}*/

/* -------------------------------------------------------- */

/* c) -------- f (x) = sqrt(x) − e^(−x) ; aproxInicial = 0.5 -----*/
/* double f(double x) {
    return sqrt(x) − pow(e,-x)
}

double df(double x) {
    return 1/(2 * sqrt(x)) + pow(e,-x)
}*/

/* -------------------------------------------------------- */

/* d) -------- f (x) = sinh(x) − sin(x) ; aproxInicial = 0 -----*/
/* double f(double x) {
    return sinh(x) - sin(x)
}

double df(double x) {
    return cosh(x) - cos(x)
}*/
