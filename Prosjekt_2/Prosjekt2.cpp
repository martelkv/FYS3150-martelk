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

using namespace std;

ofstream ofile;


void rotate(double ** A, double ** R, int k, int l, int n  );
double maxoffdiag(double ** A, int * k, int * l, int n );
mat make_matrix(int n);
void jacobi_method(double ** A, double ** R, int n);


void jacobi_method(double ** A, double ** R, int n)
{




  for ( int i = 0; i < n; i++ ) {
    for ( int j = 0; j < n; j++ ) {
  if ( i == j ) { R[i][j] = 1.0;
  }
  else { R[i][j] = 0.0;
  }
}
  }


  int k, l;
  double epsilon = 1.0e-5;
  double max_number_iterations = (double) n * (double) n * (double) n;

  int iterations = 0;
  // k er kolonne, l er rad
  double max_offdiag = maxoffdiag ( A, &k, &l, n );

  //Løkke som kaller på rotasjonsfunskjonen som roterer matrisen, og teller antall ganger dette blir gjort.
  //Fortsetter så lenge den største verdien på diagonalen er mindre enn epilon
  //Setter i og j fra maxofdiagnal til å være l og k.
  while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
    max_offdiag = maxoffdiag ( A, &k, &l, n );
    rotate (A, R, k, l, n );
    iterations++;


  }

  std::cout << "Number of iterations: " << iterations << "\n";

  return;
}


  //Maxofdiag er en funskjon som itererer gjennom elementene i matrise A. og setter diise indeksen i og j til k og l.
  // Og returner max-verdien

double maxoffdiag ( double ** A, int * k, int * l, int n )
  {
  double max = 0.0;

  for ( int i = 0; i < n; i++ ) {
    for ( int j = i + 1; j < n; j++ ) {
      if ( fabs(A[i][j]) > max ) {
        max = fabs(A[i][j]);
        //setter l og k til å være indekse i matrisen til elementet som er større enn null
        *l = i;
        *k = j;
      }
    }
    }
return max;
}


// Funskjon som roterer matrisen og finner verdiene til sin og cos
void rotate ( double ** A, double ** R, int k, int l, int n )
  {
  double s, c;
  // Rad og kolonne er byttet
  if ( A[k][l] != 0.0 ) {

    double t, tau;

    // Rad og kolonne er byttet
    tau = (A[l][l] - A[k][k])/(2*A[k][l]);
    if ( tau >= 0 ) {

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
  a_kk = A[k][k];
  a_ll = A[l][l];
// bytter matriseelementene med indekser k og l
  A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
  A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
//Hard-koder inn null-verdiene
  A[k][l] = 0.0;
  A[l][k] = 0.0;
  // Bytter om på det gjenstående elementene.
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A[i][k];
      a_il = A[i][l];
      A[i][k] = c*a_ik - s*a_il;
      A[k][i] = A[i][k];
      A[i][l] = c*a_il + s*a_ik;
      A[l][i] = A[i][l];

    }

    // Lager de nye egenvektorene
    r_ik = R[i][k];
    r_il = R[i][l];
    R[i][k] = c*r_ik - s*r_il;
    R[i][l] = c*r_il + s*r_ik;
  }

}


//Lager matrisen A
mat make_matrix (int n){

  vec analytisk = vec(n);

  vec a = vec(n-1);
  vec d = vec(n);
  vec ro = vec(n);
  ro[0] = 0;
  ro[n-1] = 1.0;

  cout <<"Dette er n:" << n<< "\n";

  double h = (ro[n-1]-ro[0])*(1./(n+1));
  cout <<"Dette er h:" << h<< "\n";


  for (int i = 0; i < n-1; i++){
    a[i] = - 1.0/(h*h);

  }
  for (int i = 0; i < n; i++){
    d[i] = 2.0/(h*h);
    ro[i] = ro[0] + (i+1)*h;
    //Regner ut den analytiske verdien
    analytisk[i] = d[i]+ 2*a[i]*cos(((i+1)* datum::pi)/(n+1));
  }
  std::cout << "Analytiske egenverdier"<< analytisk << '\n';

  //Dette er matrisen vi ønsker å rotere
  mat A  = diagmat( a, -1);
  A.diag(1) = a;
  A.diag() = d;

return A;
}



// Initialiserer matrisen
double** initMatrix(int n){

  double** B = new double*[n];

  for(int i = 0; i < n; i++){
    B[i] = new double[n];
  }

  return B;
}



//Gjør om armadillo til double **
double** matrix(mat M, int n){

  double** G = initMatrix(n);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){

      G[i][j] = M(i,j);
    }
  }
return G;
}

//Gjør om double til armadillo
mat double_to_armendillo(double** G, int n){

  mat M = mat(n,n);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){

      M(i,j) = G[i][j];
    }
  }
return M;

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
  string outfilename;

  if( argc <= 1 ){
  cout << "Bad Usage: " << argv[0] <<
    " read also output file, number of integration points and the final x values  on same line, four variables in total" << std::endl;
  exit(1);
  }

  number = atoi(argv[1]);
  outfilename=argv[2];

  int n = number;

  double** R = initMatrix(n);

  mat  midlertidig = make_matrix (n);
  //midlertidig.print("A START: ");

  double** A = matrix(midlertidig,n);



  jacobi_method(A,  R,  n);
  std::cout << "_________" << '\n';


  mat mat_eigvec = double_to_armendillo(R,n);
  mat mat_eigval = double_to_armendillo(A,n);
  vec vec_eigval = mat_eigval.diag();
  mat_eigval.print("A: ");
  mat_eigvec.print("R: ");

  /

  int indeks = vec_eigval.index_min();
  std::cout << "indeks: "<<indeks << '\n';
  std::cout << "min_verdi "<<vec_eigval(indeks) << '\n';


  vec vektor = vec(n);
  for ( int rad = 0; rad < n; rad++ ){
    vektor(rad) = mat_eigvec(rad,indeks);
  }





  //Ønsker å skrive resultatene til fil
ofile.open(outfilename);
ofile<< setprecision(8) << vektor << endl;
ofile.close();


  std::clock_t c_end = std::clock();

  double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms /1000.0 << " s\n";



std::cout  << '\n';








//TEST 1: Tester om funksjonen maxoffdiag returnerer riktig element
std::cout << "____________T_E_S_S_T________________" << '\n'<<"\n";
std::cout << "TEST1" << '\n'<<"\n";
int col;
int row;


mat testmatrise = mat (3,3);

testmatrise(0,0) = 2;
testmatrise(0,1) = 1;
testmatrise(0,2) = 10;
testmatrise(1,0) = 5;
testmatrise(1,1) = 7;
testmatrise(1,2) = 2;
testmatrise(2,0) = 1;
testmatrise(2,1) = 4;
testmatrise(2,2) = 3;
double** testmatrise1 = matrix(testmatrise,3);
std::cout << "testmatrise" << std::endl;
std::cout << testmatrise <<'\n'<<"\n";


if (maxoffdiag(testmatrise1, &col, &row, 3) == 10){
  std::cout << "col: " << col << " row: " << row << std::endl;
  std::cout << testmatrise1[row][col] << std::endl;
  std::cout << "The test passed"<<'\n';


}else{
    std::cout << "The test did not pass"<<'\n';
    std::cout << "Dette er verdien: " << maxoffdiag(testmatrise1, &col, &row, 3)<<'\n'<<"\n";

}
std::cout << "TEST 2" << '\n'<<"\n";

//TEST 2 :Tester om jacbi.metoden returnerer riktig egenverdi
mat matrise = mat(2,2);
matrise(0,0) = 3;
matrise(0,1) = 2;
matrise(1,0) = 0;
matrise(1,1) = 1;

//Har egenvektorene 1 og 3
//skjekk om dette stemmer
//Egenvekor (1,0) og (0,1)
double** A_test = matrix(matrise,2);
double** R_test = initMatrix(2);
jacobi_method(A_test,  R_test,  2);

mat mat_eigvec_test = double_to_armendillo(R_test,2);
mat mat_eigval_test = double_to_armendillo(A_test,2);

vec vec_eigval_test = mat_eigval_test.diag();



if (vec_eigval_test(0) == 0.0){
  std::cout << "The test passes" << '\n';
  }else if(vec_eigval_test(0) ==3.0){
  std::cout << "The test passes" << '\n';
    }else{
      std::cout << "Testen feiler" << '\n';

}


std::cout << "_______S_L_U_T_T___T_E_S_T_________" << '\n';



/*
  //TESTER OM MAN FÅR SAMME MATRISE TILBAKE
  mat test = make_matrix(n);
  test.print("FOR");
  double** test1 = matrix(test,n);
  mat test2 = double_to_armendillo(test1,n);
  test2.print("ETTER");

  //free_memory(R,n);
  //free_memory(A,n);




  */


return 0;
}
