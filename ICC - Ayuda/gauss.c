
/*

 * Assignatura: Introduccio a la computacio cinetifica
 * Practica 1 - Exercici 2
 * Aquest exercici resol un sistema d'equacions lineal utilitzant una matriu i resolent-la per el metode de gauss sense pivotatge de files
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int gauss(int, double **, double *, double);/* funcio que ens calcula les valors dun sistema lineal Ax = b utilitzant el metode de gauss, sense pivotatge */
void Introdueix(void); /* funcio que serveix per introduir les dades i no sobrecarregar el main */
void printvector(int, double*); /* funcio que donat un vector i la seva longitud imprimeix el vector */
void printMatrix(int,int, double **);/* funcio que donat el tamany de la matriu i la matriu, imprimeix la matriu sencera */
int comprovatoldiagonal(int ,double, double **); /* funcio que ens comprova que no hi hagin zeros a la diagonal */
double* productMatrixVector(int, double **, double *); /* funcio que calcula el producte de una matriu per un vector, retorna el punter al nou vector calculat */
double* restaVectors(int, double *, double *); /* funcio que resta dos vectors i retorna el punter del nou vector */
double modulVector(int, double *);/* funcio que calcula el modul d'un vector, com a parametres es necessita la longid del vector i el vector */
int resoltrisup(int, double **, double *, double *, double); /* funcio principal que calcula la solucions de la matriu triangular superior */

int main(void){/* funcio principal, nomes crida a la funcio Introduexi per introduir els valors necessaris */
  Introdueix();
  return 0;
}

int gauss(int n, double **A, double *v, double tol){
  int i, j,retu, pos;
  double *x, *b, *bcalculat;
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
	for ( j=i+1; j<n;++j){ 
		m = A1[j][i]/A1[i][i]; /* calculem m */
		for ( pos =0;pos <= n; ++pos){ /* multipliquem la fila anterior per m i la restem a la fila actual */
			/*if ( pos < j){ A1[j][pos] = 0.00000000;} es pot fer que a les posicions on saps que hi han d'anar zeros, posar els zeros directament i no calcular-los, però he vist que aixo feia augmentar l'error, es per aixo, que he decidit deixaro comentat.
			else{ */
			A1[j][pos] = A1[j][pos]-(m*A1[i][pos]);
			/*}*/
		}
		retu = comprovatoldiagonal(n,tol,A1); /* es comprova que la diagonal sigui diferent de zero */
		if( retu == 1){ /*si a la diagonal hi han zeros, llavors ja no podem continuar */
			printf("\n m %f \n",m);
			printf("S'ha trobat un zero a la diagonal \n "); /* enviem un  missatge derror*/ 
			return 1;
		}
		
		printf(" m %f \n",m);/* informacio de control */
		printMatrix(n,n+1,A1);
		printf("******************************* \n");
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
 
  resoltrisup(n,A2,b,x,tol); /* calculem els valors de x utilitzant la matriu triangular superior */
  bcalculat = productMatrixVector(n,A,x); /* creem el producte entre el vector x i la matriu A */
  bcalculat =restaVectors(n, bcalculat,v); /* calculem la resta entre el b calculat i x */
  modulVector(n,bcalculat); /* calculem el modul del vector resta */
  
  
  return 0;
}

int resoltrisup(int n,double **A, double *b, double *x, double tol){/* funcio principal que calcula la solucions de la matriu triangular superior */
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
  
  printf("Vector x: ");
  printvector(n,x);
  
  return 0;
}

int comprovatoldiagonal(int n, double tol, double **A){ /* funcio que ens comprova que no hi hagin zeros a la diagonal */
  int j;
  for ( j=0; j <n; ++j){/* per a cada posicio de la diagonal de la matriu */
    if (fabs(A[j][j]) <= tol){ /* mirem que el valor de la diagonal de la matriu sigui mes gran que la tolerancia */
	    printMatrix(n, n+1, A);
	    printf("\n  posicio %d %d, amb valor %f , mes petit que tolerancia %f \n", j,j,A[j][j], tol);
	    return 1; /* retornem un 1 indicant que hem trobat zeros a la diagonal */
	}
    }
  return 0; /* si n hi han zeros, retornem un 0 */
}  

void Introdueix(void){/* funcio que serveix per introduir les dades i no sobrecarregar el main */
  double **A, *b, tol; /* es declaren les variables necessaries */
  int n,i,j,scan;
 
  
  printf(" Introdueix la dimensio de la matriu(n):"); /* introduim la dimensio de la matriu per saber quants valors hem d'escoltar */
  scan = scanf("%d",&n);
  /*printf(" %d ",n);*/
  
  printf(" Introdueix la tolerancia(tol):");/* introduim el valor de la tolerancia desitjada */
  scan = scanf("%lf",&tol);
  /*printf(" %f ",tol);*/
  
  printf(" Introdueix la matriu A:"); /* introduim els valors de la matriu A, ja que ara sabem el tamany de la matriu */
  
  A = (double **)malloc(n*sizeof (double *)); /* declarem la matriu com un vector de vectors */
  if ( A == NULL){ exit(1);} /* si no hi ha memoria per declarar el vector, surtim del programa */
  
  for( i = 0; i < n; ++i){ /* per a cada columna */
    A[i] = (double *)malloc(n*sizeof(double));/* posem un vector dintre de cada columna */
    if(A[i] == NULL){ exit(2);}/* si no hi ha memoria per declarar el vector, surtim del programa */
 
    for(j = 0; j <n; ++j){
      scan = scanf("%lf",&A[i][j]); /* entrem el valor de la posicio */
      if (fabs(A[i][j]) <= tol){ A[i][j] = 0.000; }/* si el valor que introduim es menor que la tolerancia, posem un zero directament */
    }}
  
  printf(" Introdueix el vector b:"); /* introduim el vector b */
  b = (double *)malloc(n*sizeof(double));
  if ( b == NULL){ exit(1);} /* si no hi ha memoria per declarar el vector, surtim del programa */
  
  for ( i = 0; i< n; ++i){ /*per a cada posicio del vector, introduim els valors */
    scan = scanf("%lf",&b[i]);
    if ( fabs(b[i]) <= tol){ b[i] = 0.0;}/* si el valor que introduim es menor que la tolerancia, posem un zero directament */
  }
  
 /* printvector(n,b);
  printMatrix(n,A); */
  
  scan += 1; /* Ho poso per que el compilador no tregui un warning dient que la variable scan no s'utilitza*/
 
  
  printf("Matriu A\n ");
  printMatrix(n,n,A);
  printf("Vector b \n");
  printvector(n,b);
  printf("******************************* \n");
  
  gauss(n,A,b,tol);
  
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
  
  printf("Vector b calculat \n"); /* imprimim per comprovar el valor i retornem la solucio */
  printvector(n,sol);
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
  printf("Modul del vector: %f \n",solu); /* imprimim i retornem el valor del modul */
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
  