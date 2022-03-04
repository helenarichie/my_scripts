#ifndef DUST_CUDA_H
#define DUST_CUDA_H

#include<math.h>

void dust_kernel(double *dev_conserved, int nx, int ny, int nz, 
int n_ghost, int n_fields, double dt, double gamma, double *dt_array);

// general purpose functions:
void Get_Indices(int n_ghost, int nx, int ny, int nz, int &is, int &ie, int &js, int &je, int &ks, int &ke);

void Get_GTID(int &id, int &xid, int &yid, int &zid, int &tid, int nx, int ny, int nz);

double Calc_Pressure(double E, double d_gas, double vx, double vy, double vz, 
double gamma);

double Calc_Temp(double p, double n);

#ifdef DE
double Calc_Temp_DE(double d_gas, double ge, double gamma, double n);
#endif // DE

class Dust {

  public:
    double T, n, dt, d_gas, d_dust;
    double tau_sp;
    Dust(double T_in, double n_in, double dt_in, double d_gas_in, double d_dust_in) {
      T = T_in;
      n = n_in;
      dt = dt_in;
      d_gas = d_gas_in;
      d_dust = d_dust_in;
    }
    void calc_tau_sp();
    double calc_dd_dt();

  private:
    double MP = 1.6726*pow(10,-24); // proton mass in g
    double YR_IN_S = 3.154*pow(10,7); // one year in s

};

#endif // DUST_CUDA_H