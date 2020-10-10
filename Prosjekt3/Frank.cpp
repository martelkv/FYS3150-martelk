#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <math.h>
#include <stdio.h>
using namespace  std;
using namespace arma;


#include <string>
#include <vector>
//#include <utility> // std::pair

//ofstream ofile;


int n = 2; //Antall planeter
int NumberofSteps = 10000;
double FinalTime = 100;
double h = FinalTime/((double) NumberofSteps);

mat x = mat(NumberofSteps,n);
mat y = mat(NumberofSteps,n);
mat z = mat(NumberofSteps,n);
mat ax = mat(NumberofSteps,n);
mat ay = mat(NumberofSteps,n);
mat az = mat(NumberofSteps,n);
mat vx = mat(NumberofSteps,n);
mat vy = mat(NumberofSteps,n);
mat vz = mat(NumberofSteps,n);


vec m =vec(n); //Masse-veltor


double pi = acos(-1.0);
double g = 4*pi*pi;
int planetnummer = 0;

void output( mat&, mat&, mat&, string, string);



void addplanet(double x_, double y_, double z_,double vx_, double vy_, double vz_, double m_, int &planetnummer, mat &x,mat &y,mat &z,mat &vx,mat &vy, mat &vz, vec &m){

  x(0,planetnummer) = x_;    //Første rad, planetnummer er kolonne
  y(0,planetnummer) = y_;
  z(0,planetnummer) = z_;
  vx(0,planetnummer) = vx_;
  vy(0,planetnummer) = vy_;
  vz(0,planetnummer) = vz_;
  m[planetnummer] = m_;
  planetnummer++;


}

double avstand(int i, int j, int steg){ //i er planeten vi ser på, j er planeten vi blir påvirket av

  return sqrt(pow((x(steg,i)-x(steg,j)),2)+pow((y(steg,i)-y(steg,j)),2)+ pow((z(steg,i)-z(steg,j)),2));

}



void gravity(int steg){

  for(int i =0; i < planetnummer; i++){ //Planet vi ser på
    double ax_=0;
    double ay_=0;
    double az_=0;

    for(int j =0; j <planetnummer; j++){ //Hver planet rundt
      if (i != j){
        //std::cout << "Jeg skal ikke komme opp" << '\n';
        /*
        if (i ==0){
        std::cout << "Gravity " << i<<"  "<<j<< '\n';
        std::cout << x(steg,j)<< "  "<< x(steg,i)<<"    "<< avstand(i,j,steg) <<"  " << m[j]<< "   "<< steg <<"\n";
        }
        */
        ax_ += g*(x(steg,j)-x(steg,i))/pow(avstand(i,j,steg),3) *m[j];
        ay_ += g*(y(steg,j)-y(steg,i))/pow(avstand(i,j,steg),3) *m[j];
        az_ += g*(z(steg,j)-z(steg,i))/pow(avstand(i,j,steg),3) *m[j];

      }


    }
    ax(steg,i) = ax_; //akselerasjonen til planeneten med planetene rundt
    ay(steg,i) = ay_;
    az(steg,i) = az_;
    //std::cout << ax(steg,0) << '\n';
    //if (i ==0){
    //std::cout << "Gravity " << '\n';
    //std::cout << i<<"   "<<steg<< "   " << ax(steg,i) << '\n';
  //}
  }

}

void Euler(){

  for (int steg = 0; steg < NumberofSteps-1; steg++ ){
    for(int i = 0; i < planetnummer; i++ ){ //i vil si hvilekn planet vi ser på



    //r[steg][i] = sqrt(x[steg][i]*x[i]+y[i]*y[i]);

    gravity(steg);
    vx(steg+1,i) = vx(steg,i) + h*ax(steg,i);
    vy(steg+1,i) =  vy(steg,i) + h*ay(steg,i);
    vz(steg+1,i) =  vz(steg,i) + h*az(steg,i);

    x(steg+1,i) = x(steg,i) + h*vx(steg+1,i);
    y(steg+1,i) = y(steg,i) + h*vy(steg+1,i);
    z(steg+1,i) = z(steg,i) + h*vz(steg+1,i);

    }

}
//std::cout << "hei "<<x << '\n';

  //output(efile, x, y,z);   // write to file
}

void Veilet(){
  gravity(0);
  for (int steg = 0; steg < NumberofSteps-1; steg++ ){
  for(int i = 0; i < planetnummer; i++ ){ //I vil si hvilken planet vi ser på



    x(steg+1,i) = x(steg,i) + h*vx(steg,i) + (h*h/2)*ax(steg,i);
    y(steg+1,i) = y(steg,i) + h*vy(steg,i) + (h*h/2)*ay(steg,i);
    z(steg+1,i) = z(steg,i) + h*vz(steg,i) + (h*h/2)*az(steg,i);

 }


    gravity((steg+1)); //Kaller på gravity slik at den oppdaterer seg for neste seg. Vi har da regnet ut Den nye x,y og z,verdien den kommer til å bruke

    //if (i ==0 && steg < 4){
    //std::cout << x(steg,i)<< "     " << ax(steg,i)<<"    "<< ax(steg+1,i) <<"     "<<vx(steg,i)<< "      "<< steg <<"\n";
    //}
    for(int i = 0; i < planetnummer; i++ ){ //I vil si hvilken planet vi ser på

    vx(steg+1,i) = vx(steg,i) + h/2*(ax(steg+1,i)+ax(steg,i));
    vy(steg+1,i) = vy(steg,i) + h/2*(ay(steg+1,i)+ay(steg,i));
    vz(steg+1,i) = vz(steg,i) + h/2*(az(steg+1,i)+az(steg,i));

}
}
}




int main(int argc, char *argv[]) {

  //  declarations of variables
    string outfilename;
    string outfilename2;
    //int number;
    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 ){
      cout << "Bad Usage: " << argv[0] <<
        " read also output file on same line" << endl;
      //    exit(1);
    }
    else{
      outfilename=argv[1];
      outfilename2=argv[2];
    }




  //n = 2;
  addplanet(0,0,0,0,0,0,1,planetnummer, x,y,z,vx,vy,vz,m); //sola
  addplanet(1,0,0,0,2*pi,0,3e-6, planetnummer,x,y,z,vx,vy,vz,m); //earth
  //addplanet(9.641327723118710E-01,2.465760952329768E-01,7.355032830476560E-05,-4.414756238829297E-03*365,1.662854248250772E-02*365 ,-1.141917171095722E-06*365,3e-6, planetnummer,x,y,z,vx,vy,vz,m); //earth
/*
  //korigerer for massesenteret
  double x_verdi = 0;
  double y_verdi = 0;
  double z_verdi = 0;
  double vx_verdi = 0;
  double vy_verdi = 0;
  double vz_verdi = 0;
  double M = 0;


  for(int i = 0; i<planetnummer;i++){
    x_verdi += x(0,i)*m[i];
    y_verdi += y(0,i)*m[i];
    z_verdi += z(0,i)*m[i];
    x_verdi += vx(0,i)*m[i];
    y_verdi += vy(0,i)*m[i];
    z_verdi += vz(0,i)*m[i];
    M += m(i);

  }
 x_verdi /= M;
 y_verdi /= M;
 z_verdi /= M;
 vx_verdi /= M;
 vy_verdi /= M;
 vz_verdi /= M;

 for(int i = 0; i<planetnummer;i++){
   x(0,i) -= x_verdi;
   y(0,i) -= y_verdi;
   z(0,i) -= z_verdi;
   vx(0,i) -= vx_verdi;
   vy(0,i) -= vy_verdi;
   vz(0,i) -= vz_verdi;
 }

*/

  //gravity(0);


  //Euler();
  Veilet();
  output(x, y,z, outfilename, outfilename2);   // write to fil
  //std::cout << "avstand" <<avstand(0,1)<< '\n';
  //std::cout << "ax " <<ax[1]<< '\n';
  //std::cout << "ay " <<ay[1]<< '\n';
  //std::cout << "az " <<az[1]<< '\n';





  return 0;
}

void output(mat &x, mat &y, mat &z, string outfilename,string outfilename2){
ofstream file;
string outname1;
string outname2;
std::cout << "Planetnummer: "<<planetnummer << '\n';
  for(int i = 0;i<planetnummer;i++){
    std::cout << i << '\n';
    outname1= outfilename+"_"+to_string(i)+".txt";
    std::cout << outname1 << '\n';
    //outname2= outfilename2+"_"+to_string(planetnummer)+".txt";
    file.open(outname1);
    //vfile.open(outname2);

    for(int j = 0; j < NumberofSteps; j++){
      file << setiosflags(ios::showpoint | ios::uppercase);
      file << setw(15) << setprecision(8) << x(j,i);
      file<< setw(15) << setprecision(8) << y(j,i);
      file<< setw(15) << setprecision(8) << z(j,i)<<endl;
}
file.close();
}

  // close output file
//vfile.close();  // close output file
}  // end of function output
