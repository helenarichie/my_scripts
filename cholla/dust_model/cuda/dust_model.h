#ifndef DUST_CUDA_H
#define DUST_CUDA_H

#include<math.h>

typedef double Real;

void dust_update(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, int n_fields, Real dt, Real gamma, Real *dt_array);

void dust_kernel(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, int n_fields, Real dt, Real gamma, Real *dt_array);

class Dust {

  public:
    Real T_, n_, dt_, d_gas_, d_dust_;
    Real tau_sp_;
    Dust(Real T, Real n, Real dt, Real d_gas, Real d_dust) {
      T_ = T;
      n_ = n;
      dt_ = dt;
      d_gas_ = d_gas;
      d_dust_ = d_dust;
    }
    void set_tau_sp();
    Real calc_dd_dt();

  private:
    Real YR_IN_S_ = 3.154*pow(10,7); // one year in s

};

#endif // DUST_CUDA_H