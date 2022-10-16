#include <iostream>
#include <cmath>

#define LAMBDA 0.0000001

double F_even(double epsilon, double V_0, double a){
    double k = sqrt(2*(epsilon + V_0));
    double kappa = sqrt(-2*epsilon);
    return(sin(k*a/2.0)-(kappa/k)*cos(k*a/2.0));
}

double F_odd(double epsilon, double V_0, double a){
    double k = sqrt(2*(epsilon + V_0));
    double kappa = sqrt(-2*epsilon);
    return(sin(k*a/2.0)+(k/kappa)*cos(k*a/2.0));
}

double x_intercept(double V_0,double a){
    double x_0 = 0;
    double x_left = -a;
    double x_right = a;
    while(F_even(x_0)>LAMBDA){
        if(F_even(x_left)*F_even(x_0)<0){
            x_right = x_0;
        }
        else if(F_even(x_right)*F_even(x_0)<0){
            x_left = x_0;
        }
        else return x_0;

    }
    
    

}

int main(){
    return 0;
}