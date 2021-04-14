/*

 * Assignatura: Introduccio a la computacio cinetifica
 * Practica 2 - Exercici 1 part 1
 * Programa que troba zeros no lineals f(x) = 0 utilitzan el metode de newton amb una sola variable
 */

#define PI 3.141592653589793
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void Introdueix(void);
int Newton( double, double *, double, int, int);
double Fun(double , int);
double dFun(double, int);

int main(void){
  Introdueix();
  return 0;
	
}

void Introdueix(void){ /* funcio que introduexi les dades */
  double tol,x0;
  int iter, funcio, scan; 
  double *sol;
  
  sol = (double *)malloc(1*sizeof(double));
  
  printf( " 1.-f(x) = x² + sin(x)-PI \n 2.-f(x) = 1-ln(x) \n 3.-f(x)= sqrt(x)-e^(⁻x) \n 4.-f(x)=(e^x - e^(-x))/2 -sin(x) \n");/* simprimeix un menu i sintrodueixen les dades necessaries */
  funcio = 100;
  while( funcio > 4){
	printf(" Introdueix la funcio desitjada:"); 
	scan = scanf("%d",&funcio);
  }
  
  printf(" Introdueix la tolerancia(tol):");/* introduim el valor de la tolerancia desitjada */
  scan = scanf("%lf",&tol);
  
  printf(" Introdueix el nombre maxim d'iteracions:"); 
  scan = scanf("%d",&iter);
  
  printf(" Introdueix l'aproximacio inicial:");/* introduim el valor de laproximacio inicial desitjada */
  scan = scanf("%lf",&x0);

  scan = Newton(x0,sol,tol,iter,funcio); /* es crida a la funcio Newton */

  if (scan ==1){
	  printf("S'ha trobat una aproximacio prou bona \n solucio = %f \n",sol[0]);
  }
  else if ( scan == 0){
	  printf(" No s'ha trobat una bona solucio,l'ultim numero es: %f \n",sol[0]);
  }
}

int Newton(double x, double *sol, double tol, int iter, int funcio){
	
	
	while (iter > 0){/* per a cada iteracio , es calcula el seguent x tal que f(x) = 0 utilitzant el metode Newton */
		
		sol[0] = (x- (Fun(x,funcio)/dFun(x,funcio)));
		
		if ( (fabs(sol[0]-x) < tol) || (fabs(Fun(x,funcio)) < tol) ) { /* si es compleixen alguna de les seguents opcions es surt del programa */
			return 0;
		}
		if (fabs(dFun(x,funcio)) < tol){
			return 1;
		}
		x = sol[0];
		--iter;
		
	}
	return 1;
}


double Fun(double x, int funcio){ /* funcio que retorna el f(x) segona la funcio i el x */
	switch (funcio){
		case 1: 
			return pow(x,2.0)+sin(x)-PI;
			break;
		case 2:
			return 1-log(x); 
			break;
		case 3:
			return sqrt(x)-exp((-1)*x);
			break;
		case 4:
			return (exp(x)-exp((-1)*x))/2-sin(x);
			break;
		default:
			printf( " opció incorrecte");
			break;
	}
	return -1.0;
  
}
double dFun(double x, int funcio){/* funcio que retorna la derivada de  f(x) segona la funcio i el x */
	switch (funcio){
		case 1: 
			return (2*x+sin(x)); /* f(x) = 0 per x=-2.01148 i x=1.46532 */
			break;
		case 2:
			return -1/x; 
			break;
		case 3:
			return exp((-1)*x)+1/(2*sqrt(x));
			break;
		case 4:
			return (exp(-1*x)+exp(x)-2*cos(x))/2;
			break;
		default:
			printf( " opció incorrecte");
			break;
	
	}
	return -1.0;
}
