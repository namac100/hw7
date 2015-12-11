#include <cmath>
#include <iostream>

using namespace std;

void dgl(double* f, double* R, double mu);
void rk3(double* R, double mu, double dt);

int main(){
  double dt = 0.01;
  double tstart = 0.0;
  double tend = 30.0;
  double mu = 0.012277471;
  double R[4];
  
  R[0] = 0.994;
  R[1] = 0.0;
  R[2] = 0.0;
  R[3] = -2.00158510637908;
  
  for(double t = tstart; t < tend; t+=dt){
    rk3(R, mu, dt);
    cout << t << "\t" << R[0] << "\t" << R[1] << endl;
  }
  
  return 0;
}

void dgl(double* f, double* R, double mu){
  
  //Temporaere Variablen
  double r,s, r0, r1, r2, r3;

  R[0] = r0;
  R[1] = r1;
  R[2] = r2;
  R[3] = r3;
  r = sqrt(pow((r0 + mu),2.0) + pow(r1,2.0));
  s = sqrt(pow((r0 - 1 + mu),2.0) + pow(r1,2.0));
  
  //Differentialgleichung
  R[0] = r2;
  R[1] = r3;
  R[2] = r0 + 2.0*r3 - ((1 - mu)*(r0 + mu))/pow(r,3.0) - mu*r0/pow(s,3.0);
  R[3] = r1 - 2.0*r2 - ((1 - mu)*r1)/pow(r,3.0) - mu*r1/pow(s,3.0);
  
}

void rk3(double* R, double mu, double dt){  
  double temp[4], K1[4], K2[4], K3[4];
  dgl(K1, R, mu);
  temp[0] = R[0]+(dt/2.0)*K1[0];
  temp[1] = R[1]+(dt/2.0)*K1[1];
  temp[2] = R[2]+(dt/2.0)*K1[2];
  temp[3] = R[3]+(dt/2.0)*K1[3];
  dgl(K2, temp, mu);
  temp[0] = R[0] - dt*K1[0] + 2*dt*K2[0];
  temp[1] = R[1] - dt*K1[1] + 2*dt*K2[1];
  temp[2] = R[2] - dt*K1[2] + 2*dt*K2[2];
  temp[3] = R[3] - dt*K1[3] + 2*dt*K2[3];
  dgl(K3, temp, mu);
  R[0] += dt*((1/6.0)*K1[0] + (2/3.0)*K2[0] + (1/6.0)*K3[0]);
  R[1] += dt*((1/6.0)*K1[1] + (2/3.0)*K2[1] + (1/6.0)*K3[1]);
  R[2] += dt*((1/6.0)*K1[2] + (2/3.0)*K2[2] + (1/6.0)*K3[2]);
  R[3] += dt*((1/6.0)*K1[3] + (2/3.0)*K2[3] + (1/6.0)*K3[3]);
}