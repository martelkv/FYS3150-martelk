#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {
  mat A = mat(10,10);

  for (int i = 0; i < 10; i++){
    A[i,j] = 1.;
  }
  return 0;
}
