


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
#include <ctime>
#include "armadillo"
using namespace arma;
#include <iostream>
#include <cmath>
#include <math.h>
//#include "jacobi.h"

// use namespace for output and input
using namespace std;

// output fil som global variabel
ofstream ofile;


mat rotate(mat A, double ** R, int k, int l, int n  );
double maxoffdiag(mat A, int * k, int * l, int n );
mat make_matrix(int n);
void jacobi_method(mat A, double ** R, int n);


void jacobi_method(mat A, double ** R, int n)
{
  std::cout << "Jef er i jacobi_method" << '\n';
//Ta input fra kommandolinjen
//int main(int argc, char *argv[]){


  //Vi ønaker å lage en whikw-løkke som finner det største elementet i matrisen
  //og roterer matrisen slik at vi for å fjerne dette elementet


  for ( int i = 0; i < n; i++ ) {
    for ( int j = 0; j < n; j++ ) {
  if ( i == j ) { R[i][j] = 1.0;
  }
  else { R[i][j] = 0.0;
  }
}
  }
  std::cout << "R" <<R<< '\n';

  int k, l;
  double epsilon = 1.0e-8;
  double max_number_iterations = (double) n * (double) n * (double) n;
  int iterations = 0;
  double max_offdiag = maxoffdiag ( A, &k, &l, n );

  //Løkke som kaller på rotasjonsfunskjonen som roterer matrisen, og teller antall ganger dette blir gjort.
  //Fortsetter så lenge den største verdien på diagonalen er mindre enn epilon
  //Setter i og j fra maxofdiagnal til å være l og k.
  while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
    max_offdiag = maxoffdiag ( A, &k, &l, n );
    
    mat hei = rotate ( A, R, k, l, n );
    iterations++;




  }


  std::cout << "Number of iterations: " << iterations << "\n";

  return;
}


  //Lag en funskjon som itererer gjennom elementene i matrise A. og setter diise indeksen i og j til k og l.
  // Og returner max-verdien

double maxoffdiag ( mat A, int * k, int * l, int n )
  {
  double max = 0.0;
  for ( int i = 0; i < n; i++ ) {
    for ( int j = i + 1; j < n; j++ ) {
      if ( fabs(A(i,j)) > max ) {
        max = fabs(A(i,j));
        //setter l og k til å være indekse i matrisen til elementet som er større enn null
        *l = i;
        *k = j;
      }
    }
    }
return max;
}


// Function to find the values of cos and sin
mat rotate ( mat A, double ** R, int k, int l, int n )
  {
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;

    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    if ( tau > 0 ) {

      //Hvorfor er t definert som dette?
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
      }
    else {
        t = -1.0/( -tau + sqrt(1.0 + tau*tau));
      }
      c = 1/sqrt(1+t*t);
      s = c*t;
  }
  else {
    c = 1.0;
    s = 0.0;
    }
    // Deretter vil vi lage ne funksjon som bruker algorimene vi har fått oppgitt

  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
// changing the matrix elements with indices k and l A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll; A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll; A[k][l] = 0.0; // hard-coding of the zeros
  A(l,k) = 0.0;
  // and then we change the remaining elements
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);

    }

    // Finally, we compute the new eigenvectors
    r_ik = R[i][k];
    r_il = R[i][l];
    R[i][k] = c*r_ik - s*r_il;
    R[i][l] = c*r_il + s*r_ik;


  }


return A;
}



//LAGER MATRISEN A
mat make_matrix (int n){

  vec analytisk = vec(n);

  vec a = vec(n-1);
  vec d = vec(n);
  vec ro = vec(n);
  ro[0] = 0;
  ro[n] = 1.0;

  cout <<"Dette er n:" << n<< "\n";
  double h = (ro[n]-ro[0])/n;
  cout <<"Dette er h:" << h<< "\n";

  //Dette er matrisen som fylles med egenvektorene


  for (int i = 0; i < n-1; i++){
    a[i] = - 1/(h*h);
  }
  for (int i = 0; i < n; i++){
    d[i] = 2/(h*h);
    ro[i] = ro[0] + i*h;
    //Regner ut den analytiske verdien
    analytisk[i] = d[i]+ 2*a[i]*cos(((i+1)* datum::pi)/(n+1));
  }

  //cout <<"Dette er vektoren:" << a<< "\n";
  //cout <<"Dette er vektoren:" << d<< "\n";
  //Dette er matrisen vi ønsker å rotere
  mat A  = diagmat( a, -1);
  A.diag(1) = a;
  A.diag() = d;

  //cout <<"Dette er matrisen:" << A<< "\n";



return A;
}



// INITIALISERER MATRISEN
double** initMatrix(int n){

  double** B = new double*[n];

  for(int i = 0; i < n; i++){
    B[i] = new double[n];
  }

  return B;
}



//GJØR OM ARMANDILLO TIL DOUBLE **
double** matrix(mat M, int n){


  double** G = initMatrix(n);


  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){

      G[i][j] = M(i,j);


    }

  }


return G;
}


//Sletter minne til matrise
void free_memory(double**R, int n){
  for(int i = 0; i < n; i++){
    delete[] R[i];
  }
  delete[] R;

}

int main(int argc, char const *argv[]) {
  std::clock_t c_start = std::clock();
  int number;

  if( argc <= 1 ){
  cout << "Bad Usage: " << argv[0] <<
    " read also output file, number of integration points and the final x values  on same line, four variables in total" << std::endl;
  exit(1);
  }

  number = atoi(argv[1]);

  int n = number;

  double** R = initMatrix(n);

  mat  A = make_matrix (n);





  jacobi_method(A,  R,  n);
  //free_memory(R,n);
  //free_memory(A,n);


  std::clock_t c_end = std::clock();

  double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms /1000.0 << " s\n";

return 0;
}
