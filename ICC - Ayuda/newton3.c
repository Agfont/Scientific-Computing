/*

 * Assignatura: Introduccio a la computacio cinetifica
 * Practica 2 - Exercici 2
 * Programa que troba zeros no lineals f(x1,x2,x3) = 0 utilitzan el metode de newton amb tres variables
 */

#define PI 3.141592653589793
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

void Introdueix(void);
int Newton( double *, double *, double, int, int);
void Fun(double* , double *, int);
void dFun(double*, double **, int);
double randInRange( double, double);
double Max(double, double);
double Min(double, double);
double* gausspivot(int, double **, double *, double);/* funcio que ens calcula les valors dun sistema lineal Ax = b utilitzant el metode de gauss, amb pivotatge, gausspivot es igial que la funcio gauss, lo unic que se li afegeix es la funcio de pivotar files */
void printvector(int, double*); /* funcio que donat un vector i la seva longitud imprimeix el vector */
void printMatrix(int,int, double **);/* funcio que donat el tamany de la matriu i la matriu, imprimeix la matriu sencera */
int comprovatoldiagonal(int ,double, double **); /* funcio que ens comprova que no hi hagin zeros a la diagonal */
double* productMatrixVector(int, double **, double *); /* funcio que calcula el producte de una matriu per un vector, retorna el punter al nou vector calculat */
double* restaVectors(int, double *, double *); /* funcio que resta dos vectors i retorna el punter del nou vector */
double modulVector(int, double *);/* funcio que calcula el modul d'un vector, com a parametres es necessita la longid del vector i el vector */
double* resoltrisup(int, double **, double *, double *, double); /* funcio principal que calcula la solucions de la matriu triangular superior */
void pivotaFiles(int,double **,int); /* funcio que pivota les files si es necessari */

int main(void){ /* funcio principal que crida a la funcio introdueix */
  Introdueix();
  return 0;
	
}
double randInRange(double mini, double maxi){ /* funcio que donat un interval calcula un numero aleatori dintre d'aquest rang */
  sleep(1); /* li hem posat un delay per que sino l'ordinador es tant rapid
que calcula el numero aleatori utilitzant el mateix segon. 
		Com be se sap, lordinador nomes crea numeros pseudoaleatoris en funcio de linstant de temps */
  srand(time(NULL));
  return mini +  ((double)rand() / (double)RAND_MAX) * (maxi - mini);; /* crea un numero aleatori dintre del rang [mini,maxi]*/
}

double Max( double x, double y){ /* retorna el maxim de dos numeros. Se que ja  existeix una funcio max dintre la llibreria math, pero trigava menys en implementar el metode que trobar com es crida */
	if (x == y){ return x;}
	else if ( x > y ){ return x;}
	else {return y;}
}

double Min(double x, double y){ /* retorna el minim de dos numeros */
	if ( x == y){ return x;}
	else if ( x < y){ return x;}
	else { return y;}
}

void Introdueix(void){ /* funcio que introdueix els valors */
  double tol;
  int iter, funcio, scan, i;
  double *sol;
  double **cub;
  double *cantonada1;
  double *cantonada2;
  
  cantonada1 =(double *)malloc(3*sizeof(double)); /* cantonades del cub */
  cantonada2 =(double *)malloc(3*sizeof(double));
  sol = (double *)malloc(3*sizeof(double));
  
  printf( "1.- x+y+z=1 \n    y+z=0  \n x²+0.75y=0 \n \n 2.-x²+y²+z²=1 \n 1/4(x-y)²+(x+y)²+z²=1 \n (x-y)²+(x+y)²+1/4z²=1 \n"); /* simprimeix un menu per que el usuari escolleixi la funcio */
  
  printf(" Introdueix la funcio desitjada:"); 
  scan = scanf("%d",&funcio);
  
  printf(" Introdueix la tolerancia(tol):");/* introduim el valor de la tolerancia desitjada */
  scan = scanf("%lf",&tol);
  
  printf(" Introdueix el nombre maxim d'iteracions:"); 
  scan = scanf("%d",&iter);
  
  printf( " Introdueix les cantonades del cub (a1,a2,a3) i (b1,b2,b3):"); 
  for( i = 0; i < 3; ++i){
	   scan = scanf("%lf",&cantonada1[i]);
  }
  for ( i = 0; i <3; ++i){
	  scan = scanf("%lf",&cantonada2[i]);
  }
  
  /* construim el cub */
  cub = (double **)malloc(100*sizeof(double*));
  for( i=0 ; i <100; ++i){ 
	  cub[i] = (double *)malloc(3*sizeof(double));
  }
  
  
  for ( i= 0 ; i <100;++i){ /* creem 100 punts dintre del cub de forma aleatoria */
	cub[i][0] = randInRange(Min(cantonada1[0],cantonada2[0]), Max(cantonada1[0],cantonada2[0]));
	cub[i][1] = randInRange(Min(cantonada1[1],cantonada2[1]), Max(cantonada1[1],cantonada2[1]));
	cub[i][2] = randInRange(Min(cantonada1[2],cantonada2[2]), Max(cantonada1[2],cantonada2[2]));
	
	printf(" punt (x,y,z): %f %f %f ",cub[i][0], cub[i][1],cub[i][2]);
	
	scan = Newton(cub[i],sol,tol,iter,funcio);/* per a cada punt cridem a la funcio newton per a calcular el f(x1,x2,x3) = 0 */
	if (scan ==1){
		printf("S'ha trobat una aproximacio prou bona  solucio = %f %f %f \n",sol[0],sol[1],sol[2] );
	}
	else if ( scan == 0){
		printf(" No s'ha trobat una bona solucio,l'ultim numero es: %f %f %f \n",sol[0],sol[1],sol[2]);
	}
  }
    
  
}

int Newton(double *x, double *sol, double tol, int iter, int funcio){ /* funcio newton que calcula els zeros de la funcio utilitzanr el metode de newton */
	double *Ax, *f, **df;
	int i;
	int n;
	n = 3; /* dimensio matriu */
	f = (double *)malloc(n*sizeof(double)); 
	df = (double **)malloc(n*sizeof(double *));
	
	for( i = 0; i < n; ++i){
		df[i] = (double *)malloc(n*sizeof(double));
	}
	
	Fun(x,f,funcio); /* es demana lavaluacio de la funcio en el punt donat*/
	dFun(x,df,funcio);/* es demana la derivada de la funcio en el punt donat */
	
	Ax = gausspivot(n,df,f,tol); /* es crida a la funcio gausspivot per a calcular Ax */ 

	while (iter > 0){ /* per a cada iteracio es calcula el seguent x */
		for ( i = 0; i < n; ++i){
			sol[i] = x[i]-Ax[i];
			x [i] = sol[i];
		}
		
		if (modulVector(n,sol) < tol){
			return 0;
		}
		--iter;
	}
	return 1;
}


void Fun(double *x, double *f, int funcio){ /* funcio que retorna el f(x1,x2,x3) segona la funcio i el x */
	switch (funcio){
		case 1: 
			f[0] = x[0]+x[1]+x[2]-1;
			f[1] = x[1]+x[2];
			f[2] = pow(x[0],2)+0.75*x[1];
			break;
		case 2:
			f[0] = pow(x[0],2)+pow(x[1],2)+pow(x[2],2)-1;
			f[1] = 0.25*pow((x[0]-x[1]),2)+pow((x[0]+x[1]),2)+pow(x[2],2)-1;
			f[2] = pow(x[0]-x[1],2)+pow((x[0]+x[1]),2)+0.25*pow(x[2],2)-1;
			break;
		
		default:
			printf( " opció incorrecte");
			break;
	}
	
  
}

void dFun(double *x, double **df, int funcio){/* funcio que retorna la derivada de  f(x1,x2,x3) segona la funcio i el x; calcula la matriu de  derivades parcials */
	switch (funcio){
		case 1: 
			df[0][0] = 1.0;
			df[0][1] = 1.0;
			df[0][2] = 1.0;
			df[1][0] = 0.0;
			df[1][1] = 1.0;
			df[1][2] = 1.0;
			df[2][0] = 2.0*x[0];
			df[2][1] = 0.75;
			df[2][2] = 0.0;
			break;
		case 2:
			df[0][0] = 2.0*x[0];
			df[0][1] = 2.0*x[1];
			df[0][2] = 2.0*x[2];
			df[1][0] = 2.5*x[0]+1.5*x[1];
			df[1][1] = 1.5*x[0]+2.5*x[1];
			df[1][2] = 2.0*x[2];
			df[2][0] = 4.0*x[0];
			df[2][1] = 4.0*x[1];
			df[2][2] = 0.5*x[2];
			break;
		
		default:
			printf( " opció incorrecte");
			break;
	
	}
	
}

double * gausspivot(int n, double **A, double *v, double tol){
  int i, j,retu, pos;
  double *x, *b;
  double **A1,**A2; /* A1 contindra la matriu A amb el vector b */
  double m;
  
  x = (double *)malloc(n*sizeof(double));
  if(x == NULL){ exit(1);}/* si no hi ha memoria per declarar el vector, surtim del programa */
  
  A1 = (double **)malloc((n+1)*sizeof(double*)); /* num columnes*/
  if(A1 == NULL){ exit(2);}/* si no hi ha memoria per declarar el vector, surtim del programa */
	  
  A2 = (double **)malloc((n)*sizeof(double*)); /* A2 s'utilitzara mes endavant, per separar la matriu del vector b */
  if(A2 == NULL){ exit(3);}/* si no hi ha memoria per declarar el vector, surtim del programa */
	  
  b = (double *)malloc(n*sizeof(double));
  if(b == NULL){ exit(4);}/* si no hi ha memoria per declarar el vector, surtim del programa */
	  
  /* es copia la matriu A a la matriu A1 */
  
  for( i = 0; i<n;++i){
    A1[i] = (double *)malloc(n*sizeof(double)); /*num files*/
    if(A1[i] == NULL){ exit(5);}/* si no hi ha memoria per declarar el vector, surtim del programa */
    
    for ( j = 0; j<n;++j){
      A1[i][j] = A[i][j]; /* copiem A a A1 */
    }
  }
  
  /* es posa el vector b seguit de la matriu A1*/
  for (j= 0; j<n; ++j){
    A1[j][n] = v[j]; /* posem el vector b seguit de la matriu, construim una matriu adjunta */
  }
  
  
 /* comencem gauss  */
  for( i = 0; i <n; ++i){ /* per cada columna de la matriu */
	pivotaFiles(n,A1,i); /* mirem si necessitem pivotar les files */
	for ( j=i+1; j<n;++j){ 
		retu = comprovatoldiagonal(n,tol,A1); /* es comprova que la diagonal sigui diferent de zero */
		if( retu == 1){ /*si a la diagonal hi han zeros, llavors ja no podem continuar */
			/*printf("\n m %f \n",m);*/
			printf("S'ha trobat un zero a la diagonal \n "); /* enviem un  missatge derror*/ 
			/*return 1;*/
		}		
		m = A1[j][i]/A1[i][i]; /* calculem m */
		for ( pos =0;pos <= n; ++pos){ /* multipliquem la fila anterior per m i la restem a la fila actual */
			/*if ( pos < j){ A1[j][pos] = 0.00000000;} es pot fer que a les posicions on saps que hi han d'nar zeros, posar els zeros directament i no calcular-los, però he vist que aixo feia augmentar l'error, es per aixo, que he decidit deixaro comentat.
			else{ */
			A1[j][pos] = A1[j][pos]-(m*A1[i][pos]);
			/*}*/
		}
		
		
		/*printf(" m %f \n",m); informacio de control */
		/*printMatrix(n,n+1,A1);
		printf("******************************* \n");*/
	}
  }
  /*printf("\n A1[n-1][n] = %f \n ", A1[n-1][n]);*/
  
  /*separem el vector b de la matriu*/
  for ( i = 0; i<n ; ++i){ b[i]= A1[i][n];}
	
  /* separem la matriu A, obtenint la matriu triangular superior */
  for( i = 0; i < n; ++i){
	  A2[i] = (double *)malloc(n*sizeof(double));
	  if(A2[i] == NULL){ exit(8);}/* si no hi ha memoria per declarar el vector, surtim del programa */
		  
	  for ( j = 0; j<n; ++j){/* separem la matriu A sense el vector b a A2 */
		  A2[i][j] = A1[i][j];
	  }	  
  }
  
 /* printf("Matriu A2 \n ");
  printMatrix(n,n,A2);
    
  printf("Vector b \n");
  printvector(n,b); */
  
 /* printvector(n,x);*/
 
  return resoltrisup(n,A2,b,x,tol); /* calculem els valors de x utilitzant la matriu triangular superior */
}

double* resoltrisup(int n,double **A, double *b, double *x, double tol){/* funcio principal que calcula la solucions de la matriu triangular superior */
  int i,k;
  double aux;

  /* cas base Xn = Bn/Ann */
   x[n-1] = b[n-1]/A[n-1][n-1];
  /* printf( " x[n-1] = %f , b[n-1] = %f, A[n-1][n-1] = %f ", x[n-1], b[n-1],A[n-1][n-1]);*/
   
  /* cas successio Xi = (Bi-E(Aij*Xj))/Aii; */
  
  for ( i = 1; i < n; i++){ /* com que la ultima fila ja esta calculada al cas base, comencem per la penultima fila */
    aux = 0.0;
    for ( k = 0; k <=i;++k){  /* per a cada columna,a la fila actual,  anem multiplicant per el vector per la posicio corresponent, sense multiplicar la part de sota de la matriu, ja que sabem que es tot zeros */
      aux = aux + x[n-1-k]*A[n-1-i][n-1-k]; /* calculem el sumatori de la multiplicacio de la fila de la matriu per la columna del vector */
    }
    x[n-1-i] = (b[n-1-i]-aux) /A[n-1-i][n-1-i]; /* guardem els valors començant per el final x[n] */
  }
  
 /* printf("Vector x: ");
  printvector(n,x);*/
  return x;
}

int comprovatoldiagonal(int n, double tol, double **A){ /* funcio que ens comprova que no hi hagin zeros a la diagonal */
  int j;
  for ( j=0; j <n; ++j){/* per a cada posicio de la diagonal de la matriu */
    if (fabs(A[j][j]) <= tol){ /* mirem que el valor de la diagonal de la matriu sigui mes gran que la tolerancia */
	    printf(" \n ");
	   /* printMatrix(n, n+1, A);*/
	    printf("\n  posicio %d %d, amb valor %f , mes petit que tolerancia %f \n", j,j,A[j][j], tol);
	    return 1; /* retornem un 1 indicant que hem trobat zeros a la diagonal */
	}
    }
  return 0; /* si n hi han zeros, retornem un 0 */
}  

void printvector(int n, double *x){ /* funcio que donat un vector i la seva longitud imprimeix el vector */
  int i;
  for ( i = 0; i<n; ++i) {printf(" %f ", x[i]);} /* recorre tot el vector donat i imprimeix el valor de cada posicio */
  printf(" \n ");
}

void printMatrix(int n,int m, double **A){/* funcio que donat el tamany de la matriu i la mtriu, imprimeix la matriu sencera */
  int i,j;
  for ( i = 0; i < n; ++i){/* es recorre cada fila i cada columna i es va imprimint el valor de cada posicio de la matriu */
    for(j = 0; j < m; ++j){
      printf(" %f",A[i][j]);
    }
    printf(" \n "); /* despres de imprimir cada fila simprimeix un salt de linea */
  }
}

double* productMatrixVector(int n, double **A, double *x){/* funcio que calcula el producte de una matriu per un vector, retorna el punter al nou vector calculat */
  int i, j;
  double aux;
  double *sol; /* declaracio de variables */
  
  sol = (double *)malloc(n*sizeof(double)); /* aassignem la memoria al vector del producte */
  
  for ( i = 0; i <n; ++i){ /* per a cada fila de la matriu calculem la multiplicacio de la fila de la matriu per el vector */
    aux = 0.0;
    for ( j= 0; j<n; ++j){ /* calculem el sumatori de la fila de la matriu per el vector */
      aux = aux + A[i][j]*x[j];
    }
    sol[i] = aux; /* i guardem el valor al vector */
  }
  
  /*printf("Vector b calculat \n");  imprimim per comprovar el valor i retornem la solucio */
 /* printvector(n,sol);*/
  return sol;
}

double modulVector(int n, double *x){/* funcio que calcula el modul d'un vector, com a parametres es necessita la longid del vector i el vector */
  int i;
  double solu; /* declaracio de variables */
  
  solu = 0.0; 
  for ( i =0 ; i < n; ++i){ /* per a cada valor del vector, l'elebem al quadrat i el sumem als altres valors calculats */
    solu = solu + pow(x[i],2);
  }
  solu = sqrt(solu); /* fem l'arrel quadrada del modul calculat */
 /* printf("Modul del vector: %f \n",solu);  imprimim i retornem el valor del modul */
  return solu;
}

double* restaVectors(int n, double *Ax, double *b){/* funcio que resta dos vectors i retorna el punter del nou vector */
  double *sol;
  int i; /* declaracio  de variables */
  
  sol = (double *)malloc(n*sizeof(double)); /* creem el vector amb les solucions */
  
  for ( i = 0; i < n; ++i){ /* restem els elements del vector un a un */
    sol[i] = Ax[i]-b[i]; /* i els guardem al vector de solucions */
  }
  
  printf("Vector resta (Ax-b):"); /* imprimim per comprovar el resultat */
  printvector(n,sol);
  return sol; /* retornem el vector amb la diferencia entre vectors */
   
}

void pivotaFiles(int n,double **A,int columna){/* funcio que pivota les files si es necessari */
  double maxi, filaaux;
  int i, posfilaacambiar;
  
  /* busquem el maxim de la columna i guardem tambe la fila a cambiar */
  maxi = fabs(A[columna][columna]); /* iniciem pensant que la Aii actual es el mes gran */
  posfilaacambiar = columna;
  
  for ( i = columna+1; i < n; ++i){ /* si trobem algun membre de la columna mes gran que el Aii,  ho guardem */
    if (fabs(A[i][columna]) > maxi){
      posfilaacambiar = i;
      maxi = fabs(A[i][columna]);
    }
  }
  
  /*filaaux = (double *)malloc((n+1)*sizeof(double));  punter auxiliar que guardara un vector mentre movem les files */
  if ( posfilaacambiar != columna){
    
    for ( i = 0; i < n+1; ++i){ /* pivotem les files */
      filaaux = A[columna][i]; /* posem la fila a "baixar" en un auxiliar */
      A[columna][i] = A[posfilaacambiar][i]; /* baixem la fila */
      A[posfilaacambiar][i] = filaaux;/* posem la fila que estava en l'auxiliar a la matriu */
     
    }
   /* printf(" Matriu pivotada: fila %d <-> fila %d \n",columna,posfilaacambiar);
    printMatrix(n,n+1,A);*/
  }
  
}