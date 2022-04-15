#ifndef DUST_CUDA_H
#define DUST_CUDA_H

#include<math.h>

typedef double Real;

void Dust_Update(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, int n_fields, Real dt, Real gamma);

__global__ void Dust_Kernel(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, int n_fields, Real dt, Real gamma);

void Conserved_Init(Real *host_conserved, Real rho, Real vx, Real vy, Real vz, Real P, Real rho_dust, Real gamma, int n_cells, int nx, int ny, int nz, int n_ghost, int n_fields);

class Dust {

  public:
    Real T_, n_, dt_, d_gas_, d_dust_;
    Real tau_sp_;
    __device__ Dust(Real T, Real n, Real dt, Real d_gas, Real d_dust) {
      T_ = T;
      n_ = n;
      dt_ = dt;
      d_gas_ = d_gas;
      d_dust_ = d_dust;
    }
    __device__ void set_tau_sp();
    __device__ Real calc_dd_dt();

  private:
    Real YR_IN_S_ = 3.154*pow(10,7); // one year in s

};

#endif // DUST_CUDA_H