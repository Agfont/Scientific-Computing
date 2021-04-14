/*

 * Assignatura: Introduccio a la computacio cinetifica
 * Practica 2 - Exercici 1 part 2
 * Programa que troba les conques datraccio de zeros en polinomis utilitzant el metode de Newton
 */

#define PI 3.141592653589793
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void Introdueix(void);
int Newton( double, double *, double, int, int);
double Fun(double , int);
double dFun(double, int);
double M(int);

int main(void){
  Introdueix();
  return 0;
}

void Introdueix(void){/* funcio que introduexi les dades */
  double tol, m, h;
  int iter, funcio, scan; 
  double *sol;
  
  sol = (double *)malloc(1*sizeof(double));
  
  printf( "1.- P(x) = x²-1 \n 2.-P(x) = x³-x \n 3.-P(x)= 3x³-x+1 \n 4.-P(x)=x⁴+1 \n 5.-P(x) = Prod(i=1,6)(x-10i/i+1)  \n 6.-P(x) = Prod(i=1,6)(x-(i+1)/10i) \n");/* simprimeix un menu i sintrodueixen les dades necessaries */
  
  printf(" Introdueix la funcio desitjada:"); 
  scan = scanf("%d",&funcio);
  
  printf(" Introdueix la tolerancia(tol):");/* introduim el valor de la tolerancia desitjada */
  scan = scanf("%lf",&tol);
  
  printf(" Introdueix el nombre maxim d'iteracions:"); 
  scan = scanf("%d",&iter);
  
  printf(" Introdueix el pas h:");/* introduim */
  scan = scanf("%lf",&h);
  
  m = M(funcio); /* es crida el metode M segons la funcio seleccionada */
  
  m = (-1.0)*m; /* comencem desde -M fins a M */
  while( m <= M(funcio) ){
	printf( " M :%f",m);
	scan = Newton(m,sol,tol,iter,funcio);
	if (scan ==1){
		printf("S'ha trobat una aproximacio prou bona \n solucio = %f \n",sol[0]);
	}
	else if ( scan == 0){
		printf(" No s'ha trobat una bona solucio,l'ultim numero es: %f \n",sol[0]);
	}
	m = m + h; /* sincrementa m segons h */
	}
	
}

int Newton(double x, double *sol, double tol, int iter, int funcio){
	
	while (iter > 0){/* per a cada iteracio , es calcula el seguent x tal que f(x) = 0 utilitzant el metode Newton */
		sol[0] = (x- (Fun(x,funcio)/dFun(x,funcio)));
		
		if ( (fabs(sol[0]-x) < tol) || (fabs(Fun(x,funcio)) < tol) ) {/* si es compleixen alguna de les seguents opcions es surt del programa */
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


double Fun(double x, int funcio){/* funcio que retorna el f(x) segona la funcio i el x */
	switch (funcio){
		case 1: 
			return (pow(x,2.0)-1.0);
			break;
		case 2:
			return (pow(x,3.0)-x); 
			break;
		case 3:
			return (3.0*(pow(x,3.0))-x+1.0);
			break;
		case 4:
			return ((pow(x,4.0))+1.0);
			break;
		case 5:
			return ((pow(x,6.0))-(617/14)*(pow(x,5.0))+(50705/63)*(pow(x,4.0))-(981775/126)*(pow(x,3.0))+(2655500/63)*(pow(x,2.0))-(845000/7)*x+(1000000/7));
			break;
		case 6:return ((pow(x,6.0))-(169/200)*(pow(x,5.0))+(5311/18000)*(pow(x,4.0))-(39271/720000)*(pow(x,3.0))+(10141/1800000)*(pow(x,2.0))+(617/2000000)*x+(7/1000000));	
			break;
		default:
			printf( " opció incorrecte");
			break;
	}
	return -1.0;
  
}
double dFun(double x, int funcio){/* funcio que retorna la derivada de  f(x) segons la funcio i el x */
	switch (funcio){
		case 1: 
			return (2.0*x);
			break;
		case 2:
			return (3.0*(pow(x,2.0))-1.0); 
			break;
		case 3:
			return (9.0*(pow(x,2.0))-1.0);
			break;
		case 4:
			return (4.0*(pow(x,3.0)));
			break;
		case 5:
			return ((756*(pow(x,5.0))- 27765*(pow(x,4.0)) + 405640*(pow(x,3.0))-2945325*(pow(x,2.0))-10622000*x-15210000)/126);
			break;
		case 6:
			return (6*(pow(x,5.0))-(169/40)*(pow(x,4.0))+(5311/4500)*(pow(x,3.0))-(39271/240000)*(pow(x,2.0))+(10141/900000)*x +(617/2000000));
			break;
		default:
			printf( " opció incorrecte");
			break;
	
	}
	return -1.0;
}
double M(int funcio){ /* retorna la M calculada */
	switch (funcio){
		case 1: 
			return 2.0;
			break;
		case 2:
			return 2.0; 
			break;
		case 3:
			return (5.0/3.0);
			break;
		case 4:
			return 2.0;
			break;
		case 5:
			return 38679.21428571428;
			break;
		case 6:
			return 1.900034240740741;
			break;
		default:
			printf( " opció incorrecte");
			break;
	
	}
	return -1.0;
}

