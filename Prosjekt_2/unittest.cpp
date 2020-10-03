
//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
//#include "catch.hpp"
#include "Prosjekt2.cpp"

int main(int argc, char const *argv[]) {

int* k;
int* l;

mat testmatrise = mat (3,3);

testmatrise(0,0) = 2;
testmatrise(0,1) = 1;
testmatrise(0,2) = 1;
testmatrise(1,0) = 5;
testmatrise(1,1) = 7;
testmatrise(1,2) = 2;
testmatrise(2,0) = 1;
testmatrise(2,1) = 4;
testmatrise(2,2) = 3;
double** testmatrise1 = matrix(testmatrise,3);
std::cout << "testmatrise" <<testmatrise <<'\n';
//TEST_CASE("Hente riktig maksverdi"){

  //REQUIRE(max_offdiag(testmatrise1,&k, &l,3)==5) // Må man ha med l og k

if (maxoffdiag(testmatrise1,k, l,3)==5){
  std::cout << "The test passed"<<'\n';

}else{
    std::cout << "The test did not pass"<<'\n';

}


  return 0;
}
//Kan også sjekke ortogonalite. Ganger man to basiser sammen skal de bli 1. Ganger man to like sammen skal det bli 0.

/*
mat matrise = mat(2,2);
matrise(0,0) = 3;
matrise(0,1) = 2;
matrise(1,0) = 0;
matrise(1,1) = 1;

//Har egenvektorene 1 og 3
//skjekk om dette stemmer
//Egenvekor (1,0) og (0,1)
double** A = matrix(matrise,n);
double** R = initMatrix(n);
jacobi_method(A,  R,  n);
mat mat_eigvec = double_to_armendillo(R,n);
mat mat_eigval = double_to_armendillo(A,n);
vec vec_eigval = mat_eigval.diag();

TEST_CASE("Gir sen riktig egenverdi"){

  REQUIRE(vec_eigval == 3 eller 1) //
}
*/
