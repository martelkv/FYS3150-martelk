#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
#include <ctime>
#include <armadillo>
#include <cmath>


using namespace std;
using namespace arma;

// output fil som global variabel
ofstream ofile;


//Dette programmet regner ut egenverdiene for en matrise ved bruk av Jacobi-metoden.
//Koden er ganske lik "Prosjekt2.cpp", men make_matrix er endret og har gjort om noen av funskjonene
// til å ta inn aramdillo mat istenfor double**. Har også fjernet matrisen R, ettersom vi her ikke er like
// interessert i å se på vektorene.


mat rotate(mat A, int k, int l, int n  );
double maxoffdiag(mat A, int * k, int * l, int n );
mat make_matrix(int n);
mat jacobi_method(mat A, int n);


mat jacobi_method(mat A, int n)
{



  int k, l;
  double epsilon = 1.0e-8;
  double max_number_iterations = 3*(double) n * (double) n;
  int iterations = 0;
  double max_offdiag = maxoffdiag ( A, &k, &l, n );

  //Løkke som kaller på rotasjonsfunskjonen som roterer matrisen, og teller antall ganger dette blir gjort.
  //Fortsetter så lenge den største verdien på diagonalen er mindre enn epilon
  //Setter i og j fra maxofdiagnal til å være l og k.
  while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
    max_offdiag = maxoffdiag ( A, &k, &l, n );

    A = rotate ( A, k, l, n );
    iterations++;


  }

  std::cout << "Number of iterations: " << iterations << "\n";

  return A;
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


// Funskjon som roterer matrisen og finner verdiene til sin og cos
mat rotate ( mat A, int k, int l, int n )
  {
  double s, c;
  // Rad og kolonne er byttet
  if ( A(k,l) != 0.0 ) {
    double t, tau;

    // Rad og kolonne er byttet
    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    if ( tau > 0 ) {
      //t er definert slik for å unngå "loss of nummerical precison"
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
  // Lager en funksjon som bruker utrykkene vi har fått oppgitt
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  // bytter matriseelementene med indekser k og l
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  //Hard-koder inn null-verdiene
  A(k,l) = 0.0;
  A(l,k) = 0.0;
  // Bytter om på det gjenstående elementene.
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);

    }

  }
return A;
}




//Lager matrisen for et elektron
mat make_matrix (int n){

  std::cout << "n"<<n << '\n';
  vec d = vec(n);
  vec e = vec(n-1);
  vec ro = vec(n);
  vec v = vec(n-1);
  ro[0] = 0;
  ro[n-1] = 5.;

  double h = (ro[n-1]-ro[0])*(1./(n+1));

  //Dette er matrisen som fylles med egenvektorene
  for (int i = 0; i < n; i++){

    ro[i] = ro[0] + (i+1)*h;
    v[i] = ro[i]*ro[i];
    d[i] = 2/(h*h) + v[i];
  }

  for (int i = 0; i < n-1; i++){
    e[i] = - 1/(h*h);

  }


  mat A  = diagmat( e, -1);

  A.diag(1) = e;

  A.diag() = d;


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

//Lager matrisen for 2 elektoner
mat make_matrix_two_electons (int n, double w){

  vec analytisk = vec(n);

  std::cout << "n"<<n << '\n';
  vec d = vec(n);
  vec a = vec(n-1);
  vec ro = vec(n);
  vec v = vec(n-1);
  ro[0] = 0;
  ro[n-1] = 5.;

  double h = (ro[n-1]-ro[0])*(1./(n+1));

  //Dette er matrisen som fylles med egenvektorene
  for (int i = 0; i < n; i++){

    ro[i] = ro[0] + (i+1)*h;
    v[i] = ro[i]*ro[i]*w*w + (1/ro[i]);
    d[i] = 2/(h*h) + v[i];

  }
  for (int i = 0; i < n-1; i++){
    a[i] = - 1.0/(h*h);

  }


  mat A  = diagmat( a, -1);
  A.diag(1) = a;
  A.diag() = d;

return A;
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


  //Regner ut egenverdiene for et elektron
  std::cout << '\n'<<"____ET_ELEKTROEN_______"<< "\n";


  mat A = make_matrix (n);
  A = jacobi_method(A,  n);


  vec vec_eigvals = vec(n);
  for (int i = 0; i < n; i++) vec_eigvals[i] = A(i,i);
  vec_eigvals = sort(vec_eigvals, "ascend");

  for (int i = 0; i< 4; i++){
    cout << vec_eigvals(i) << endl;
  }


  //Regner ut egenverdiene for 2 elektron
  std::cout << '\n'<<"____TO_ ELEKTROER_______"<< "\n";

  double w1 = 0.01;
  double w2 = 0.5;
  double w3 = 1;
  double w4 = 5;

  mat to_elektoner1 = make_matrix_two_electons(n,w1);
  mat to_elektoner2 = make_matrix_two_electons(n,w2);
  mat to_elektoner3 = make_matrix_two_electons(n,w3);
  mat to_elektoner4 = make_matrix_two_electons(n,w4);



  to_elektoner1 = jacobi_method(to_elektoner1,n);
  to_elektoner2 = jacobi_method(to_elektoner2,n);
  to_elektoner3 = jacobi_method(to_elektoner3,n);
  to_elektoner4 = jacobi_method(to_elektoner4,n);



  vec  to_elektoner1_egenverdier= vec(n);
  vec  to_elektoner2_egenverdier= vec(n);
  vec  to_elektoner3_egenverdier= vec(n);
  vec  to_elektoner4_egenverdier= vec(n);



  for (int i = 0; i < n; i++) {
    to_elektoner1_egenverdier[i] = to_elektoner1(i,i);
    to_elektoner2_egenverdier[i] = to_elektoner2(i,i);
    to_elektoner3_egenverdier[i] = to_elektoner3(i,i);
    to_elektoner4_egenverdier[i] = to_elektoner4(i,i);
  }

  to_elektoner1_egenverdier = sort(to_elektoner1_egenverdier, "ascend");
  to_elektoner2_egenverdier = sort(to_elektoner2_egenverdier, "ascend");
  to_elektoner3_egenverdier = sort(to_elektoner3_egenverdier, "ascend");
  to_elektoner4_egenverdier = sort(to_elektoner4_egenverdier, "ascend");



  for (int i = 0; i< 4; i++){
    cout << to_elektoner1_egenverdier(i)<< "  "<< to_elektoner2_egenverdier(i)<< "   "<< to_elektoner3_egenverdier(i)<< "  "<< to_elektoner4_egenverdier(i)<< "  " << endl;

}



  std::clock_t c_end = std::clock();

  double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms /1000.0 << " s\n";

return 0;
}
