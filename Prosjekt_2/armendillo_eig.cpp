
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
//#include "jacobi.h"


// use namespace for output and input
using namespace std;

// output fil som global variabel
ofstream ofile;


//LAger matrisen A som vi skal finne egenverdiene og og egenvektorene til.
mat make_matrix(int n){
  std::cout << "Jeg er i make_matrix" << '\n';
  vec a = vec(n-1);
  vec d = vec(n);
  vec ro = vec(n);
  ro[0] = 0;
  ro[n-1] = 1.0;

  cout <<"Dette er n:" << n<< "\n";

  double h = (ro[n-1]-ro[0])*(1./(n+1));
  cout <<"Dette er h:" << h<< "\n";

  //Dette er matrisen som fylles med egenvektorene

  for (int i = 0; i < n-1; i++){
    a[i] = - 1/(h*h);
  }
  for (int i = 0; i < n; i++){
    d[i] = 2/(h*h);
    ro[i] = ro[0] + (i+1)*h;
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

  mat A = make_matrix(n);



  std::cout << A.is_diagmat() << '\n';
  //Oppretter vektorene og matrisen egenverdiene og egenvektorene skal befinne seg i.
  vec eigval(n);
  mat eigvec(n,n);

  //Finner egenvektorene og egenverdiene i matrisen A
  eig_sym(eigval,eigvec, A);
  //eigval.print();
  std::cout  << '\n';
  //eigvec.print();

  //Finner minste indeks i vektoren med egenverdier
  int indeks = eigval.index_min();
  std::cout << "indeks: "<<indeks << '\n';
  std::cout << "min_verdi "<<eigval(indeks) << '\n';

  //Finner egenvektoren som tilhører den minste egenverdien
  vec vektor = vec(n);
  for ( int rad = 0; rad < n; rad++ ){
    vektor(rad) = eigvec(rad,indeks);
    }
    //Egenvektoren til den minste egenverdien
  //std::cout << vektor << '\n';

  std::clock_t c_end = std::clock();

  double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms /1000.0 << " s\n";

  return 0;
}
