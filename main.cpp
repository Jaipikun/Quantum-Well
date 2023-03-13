#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
#define V_0_min 0
#define V_0_max 10   // g³êbokoœæ studni potencja³u - energia potencjalna
#define A_min 0
#define A_max 10    // szerokoœæ studni potencja³u
#define Step_V_0 0.01
#define Step_A 0.01
#define EPSILON_E_H 0.005
//#define EPSILON V_0 / 100000.0 // precyzja wyznaczenia energii
//#define DELTA_E V_0 / 100.0  // przedzia³y energii, w ka¿dym szukamy max 1 miejsca zerowego
///////////////////////////////////////////////////////////////////////////////

fstream plik_out_1;
fstream plik_out_2;

///////////////////////////////////////////////////////////////////////////////

fstream data_out_1;
fstream data_out_2;

///////////////////////////////////////////////////////////////////////////////

double Zero(double x_min, double x_max, double epsilon, double V_0, double A , double (*Function)(double,double,double))
{
  double x_medium;

  if( Function(x_min,V_0,A) == 0 ) return x_min;
  
  if( Function(x_max,V_0,A) == 0 ) return x_max;

  while( (x_max-x_min) > epsilon)
  {
    
    x_medium = (x_max+x_min)/2.0;
    
    if( Function(x_min,V_0,A)*Function(x_medium,V_0,A) < 0 )
    {
       x_max = x_medium;
    }

    else 
    {
      x_min = x_medium;
    }
  }

  return x_medium;
}

double F_even( double E, double V_0, double A )
{
  double k = sqrt( 2.0*(V_0+E) );
  double F = sin( k*A/2.0 ) - cos( k*A/2.0 ) * sqrt( -2.0*E )/k;

  return F;
}

double F_odd( double E, double V_0, double A )
{
  double k = sqrt( 2.0*(V_0+E) );
  double F = sin( k*A/2.0 ) + cos( k*A/2.0 ) * k/sqrt( -2.0*E );

  return F;
}

void Graph(double V_0,double A,double DELTA_E)
{
  double E = -V_0 + 0.001*V_0;
  
  while( E < 0 )
  {
    data_out_1 <<E<<"\t"<<F_even(E,V_0,A)<<"\t"<< F_odd(E,V_0,A)<<"\n";
    E = E + DELTA_E;
  }

}



///////////////////////////////////////////////////////////////////////////////
int main()
{
    
    double V_0 = V_0_min; //starting value of V_0
    double A , E , E_0;

    double n = abs(V_0_max-V_0_min)/Step_V_0;   // amount of V_0 values to be checked

  //file output
    data_out_1.open("out_1.txt", ios::out);
    data_out_2.open("out_2.txt", ios::out);
    data_out_1.precision(7);      
    data_out_2.precision(10);     


    
    for(int i = 0 ; i<n;i++)
    {
        A = A_min; //starting value of A

        double DELTA_E = V_0 / 100.0;     //  energy range for x-intercept search
        double EPSILON = V_0 / 100000.0;  //  precision of energy's value

        double m = abs(A_max-A_min)/Step_A;   // amount of A values to be checked per singular value of V_0
        if(V_0!=0)
        {
        for(int j = 0;j<m;j++)
            {
            
            bool E_1 = false; // Energy level 1 checker
            bool E_2 = false; // Energy level 2 checker
            int counter = 0;  // Current energy level


            E = -V_0 + 0.001*V_0;
            
            while( E < 0 )
              {
                double E_H = -1/(2*pow(counter+1,2)); // Hydrogen level

                if( F_even(E,V_0,A)*F_even(E+DELTA_E,V_0,A) < 0 )
                {
                  counter++;
                  E_0 = Zero(E, E+DELTA_E, EPSILON,V_0,A, F_even);

<<<<<<< HEAD
                if( F_even(E,V_0,A)*F_even(E+DELTA_E,V_0,A) < 0 ){
                counter++;
                E_0 = Zero(E, E+DELTA_E, EPSILON,V_0,A, F_even);
                  if(E_0>=E_H-(EPSILON_E_H/pow(10,counter-1)) && E_0<=E_H+(EPSILON_E_H/pow(10,counter-1)))
=======
                  if(E_0>=E_H-EPSILON_E_H && E_0<=E_H+EPSILON_E_H)
                  {
                    if(E_1 || E_H==-1.0/2.0)
>>>>>>> af11f60fcdf2ca75ce5403938912d5f1cc0e0e02
                    {
                      if(!E_1)  
                      {
<<<<<<< HEAD
                        if(!E_1){
                          plik_out_2<<"\n\nV_0: "<<V_0<<" A: "<<A<<endl;
                          
                          E_1 = true;
                        }
                        else if(!E_2 && E_H==-1.0/8.0){
                          plik_out_1<<"\n\n["<<i<<":"<<j<<"]V_0:  "<<V_0<<"  A: "<<A<<"\n\n";
                          Graph(V_0,A,DELTA_E);
                          E_2 = true;
                        }
                        if(E_1 && (E_2 || E_H==-1.0/2.0)){
                          plik_out_2 << counter << " " << E_0 << " " << EPSILON_E_H/pow(10,counter-1) << endl;
                        }
                        
=======
                        data_out_2<<"\n\nV_0: "<<V_0<<" A: "<<A<<endl;    
                        E_1 = true;
>>>>>>> af11f60fcdf2ca75ce5403938912d5f1cc0e0e02
                      }
                      else if(!E_2 && E_H==-1.0/8.0)  
                      {
                        data_out_1<<"\n\n["<<i<<":"<<j<<"]V_0:  "<<V_0<<"  A: "<<A<<"\n\n";
                        Graph(V_0,A,DELTA_E);
                        E_2 = true;
                      }
                      if(E_1 && (E_2 || E_H==-1.0/2.0))
                      {
                        data_out_2 << counter << " " << E_0 << endl;
                      }
                         
                    }
                  }

                  
                }

                E_H = -1/(2*pow(counter+1,2));

<<<<<<< HEAD
                if( F_odd(E,V_0,A)*F_odd(E+DELTA_E,V_0,A) < 0 ){
                counter++;
                E_0 = Zero(E, E+DELTA_E, EPSILON,V_0,A, F_odd);
                if(E_0>=E_H-(EPSILON_E_H/pow(10,counter-1)) && E_0<=E_H+(EPSILON_E_H/pow(10,counter-1)))
=======
                if( F_odd(E,V_0,A)*F_odd(E+DELTA_E,V_0,A) < 0 )
                {
                  counter++;
                  E_0 = Zero(E, E+DELTA_E, EPSILON,V_0,A, F_odd);
                  
                  if(E_0>=E_H-EPSILON_E_H && E_0<=E_H+EPSILON_E_H)
                  {
                    if(E_1 || E_H==-1.0/2.0)
>>>>>>> af11f60fcdf2ca75ce5403938912d5f1cc0e0e02
                    {
                      if(!E_1)
                      {
<<<<<<< HEAD
                        if(!E_1){
                          plik_out_2<<"\n\nV_0: "<<V_0<<" A: "<<A<<endl;
                          E_1 = true;
                        }
                        else if(!E_2 && E_H==-1.0/8.0){
                          plik_out_1<<"\n\n["<<i<<":"<<j<<"]V_0:  "<<V_0<<"  A: "<<A<<"\n\n";
                          Graph(V_0,A,DELTA_E);
                          E_2 = true;
                        }
                        if(E_1 && (E_2 || E_H==-1.0/2.0)){
                          plik_out_2 << counter << " " << E_0 << " " << EPSILON_E_H/pow(10,counter-1) << endl;
                        }
=======
                        data_out_2<<"\n\nV_0: "<<V_0<<" A: "<<A<<endl;
                        E_1 = true;
                      }

                      else if(!E_2 && E_H==-1.0/8.0)
                      {
                        data_out_1<<"\n\n["<<i<<":"<<j<<"]V_0:  "<<V_0<<"  A: "<<A<<"\n\n";
                        Graph(V_0,A,DELTA_E);
                        E_2 = true;
                      }

                      if(E_1 && (E_2 || E_H==-1.0/2.0))
                      {
                        data_out_2 << counter << " " << E_0 << endl;
>>>>>>> af11f60fcdf2ca75ce5403938912d5f1cc0e0e02
                      }
                    }
                  }
                }

                E = E + DELTA_E;  // next E value
              }

            
            A+=Step_A;  // next A value
            }
        }
        V_0+=Step_V_0;  // next V_0 value
    }
    
    
    data_out_1.close();
    data_out_2.close();
    
    
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
