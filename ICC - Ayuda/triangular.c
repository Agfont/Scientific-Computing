
/*

 * Assignatura: Introduccio a la computacio cinetifica
 * Practica 1 - Exercici 1
 * Aquest exercici resol una matri triangular superior(zeros sota la diagonal ) i calcula el el vector diferencia entre la solucio 
 * calculada i la real
 * */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int resoltrisup(int, double **, double *, double *, double); /* funcio principal que calcula la solucions de la matriu triangular superior */
void Introdueix(void); /* funcio que serveix per introduir les dades i no sobrecarregar el main */
void printvector(int, double*); /* funcio que donat un vector i la seva longitud imprimeix el vector */
void printMatrix(int, double **);/* funcio que donat el tamany de la matriu i la mtriu, imprimeix la matriu sencera */
double* productMatrixVector(int, double **, double *); /* funcio que calcula el producte de una matriu per un vector, retorna el punter al nou vector calculat */
double* restaVectors(int, double *, double *); /* funcio que resta dos vectors i retorna el punter del nou vector */
double modulVector(int, double *);/* funcio que calcula el modul d'un vector, com a parametres es necessita la longid del vector i el vector */

int main(void){  /* funcio principal, nomes crida a la funcio Introduexi per introduir els valors necessaris */
  Introdueix();
  return 0;
}

void printvector(int n, double *x){ /* funcio que donat un vector i la seva longitud imprimeix el vector */
  int i;
  for ( i = 0; i<n; ++i) {printf(" %f ", x[i]);} /* recorre tot el vector donat i imprimeix el valor de cada posicio */
  printf(" \n ");
}

void printMatrix(int n, double **A){/* funcio que donat el tamany de la matriu i la mtriu, imprimeix la matriu sencera */
  int i,j;
  for ( i = 0; i < n; ++i){/* es recorre cada fila i cada columna i es va imprimint el valor de cada posicio de la matriu */
    for(j = 0; j<n; ++j){
      printf(" %f",A[i][j]);
    }
    printf(" \n "); /* despres de imprimir cada fila simprimeix un salt de linea */
  }
}

void Introdueix(void){/* funcio que serveix per introduir les dades i no sobrecarregar el main */
  double **A, *b,*x, tol; /* es declaren les variables necessaries */
  int n,i,j,scan;
  
  x = (double *)malloc(n*sizeof(double)); /* vector on estara el vector x amb els resultats */
  
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
  
  resoltrisup(n,A,b,x,tol); /* cridem a la funcio resoltripsup per que resolgui la matriu triangular superior, amb el tamany de la matriu, , la matriu A, el vector b, el vector de solucions x, i la tolerancia */ 
  }

int resoltrisup(int n,double **A, double *b, double *x, double tol){/* funcio principal que calcula la solucions de la matriu triangular superior */
  int i,k;
  double aux;
  double *sol; /* declaracio de variables */
  
   /* cas base Xn = Bn/Ann */
  x[n-1] = b[n-1]/A[n-1][n-1];
  
  /* cas successio Xi = (Bi-E(Aij*Xj))/Aii; */
  
  for ( i = 1; i < n; i++){ /* com que la ultima fila ja esta calculada al cas base, comencem per la penultima fila */
    aux = 0.0;
    for ( k = 0; k <=i;++k){  /* per a cada columna,a la fila actual,  anem multiplicant per el vector per la posicio corresponent, sense multiplicar la part de sota de la matriu, ja que sabem que es tot zeros */
      aux = aux + x[n-1-k]*A[n-1-i][n-1-k]; /* calculem el sumatori de la multiplicacio de la fila de la matriu per la columna del vector */
    }
    x[n-1-i] = (b[n-1-i]-aux) /A[n-1-i][n-1-i]; /* guardem els valors començant per el final x[n] */
  }
  
  printf("Vector x \n");
  printvector(n,x);
  
  /* calculem el modul de Ax-b utilitzant les funcions auxiliars */
  sol = productMatrixVector(n,A,x);
  sol = restaVectors(n,sol,b);
  aux = modulVector(n,sol);
  
  return 0;
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
  
  