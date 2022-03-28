#include"cpp_dust_kernel.h"
#include<math.h>
#include<vector>

#define TINY_NUMBER 1.0e-20
#define PI 3.141592653589793
#define MP 1.672622e-24 // mass of proton, grams
#define KB 1.380658e-16 // boltzmann constant, cgs

#define TIME_UNIT 3.15569e10 // 1 kyr in s
#define LENGTH_UNIT 3.08567758e21 // 1 kpc in cm
#define MASS_UNIT 1.98855e33 // 1 solar mass in grams
#define DENSITY_UNIT (MASS_UNIT/(LENGTH_UNIT*LENGTH_UNIT*LENGTH_UNIT))
#define VELOCITY_UNIT (LENGTH_UNIT/TIME_UNIT)
#define ENERGY_UNIT (DENSITY_UNIT*VELOCITY_UNIT*VELOCITY_UNIT)
#define PRESSURE_UNIT (DENSITY_UNIT*VELOCITY_UNIT*VELOCITY_UNIT)
#define SP_ENERGY_UNIT (VELOCITY_UNIT*VELOCITY_UNIT)

#define threadIdx 1
#define blockIdx 1
#define blockIdy 1
#define gridDimx 1
#define blockDimx 1


int main() {
    double dev_conserved[6] = {1, 2, 3, 4, 5, 6};
    int nx = 1;
    int ny = 1;
    int nz = 1;
    int n_ghost = 0;
    int n_fields = 0;
    double dt = 10000;
    double gamma = 5/2;
    double dt_array[6] = {dt, dt, dt, dt, dt, dt};

    dust_kernel(dev_conserved, nx, ny, nz, n_ghost, n_fields, dt, gamma, dt_array);
}

void dust_kernel(double *dev_conserved, int nx, int ny, int nz, int n_ghost, 
int n_fields, double dt, double gamma, double *dt_array) {
   
    // get grid indices
    int n_cells = nx * ny * nz;
    int is, ie, js, je, ks, ke;
    Get_Indices(n_ghost, nx, ny, nz, is, ie, js, je, ks, ke);

    // get a global thread ID
    int id;
    int xid, yid, zid;
    int tid;
    Get_GTID(id, xid, yid, zid, tid, nx, ny, nz);
  
    // define physics variables
    double d_gas, d_dust; // fluid mass densities
    double n; // gas number density
    double T, E, p; // temperature, energy, pressure
    double mu = 0.6; // mean molecular weight
    double vx, vy, vz; // velocities
    #ifdef DE
    double ge;
    #endif // DE

    // define integration variables
    double dd_dt; // instantaneous rate of change in dust density
    double dd; // change in dust density at current time-step
    double dd_max = 0.01; // allowable percentage of dust density increase
    double dt_sub; //refined timestep
    
    if (xid >= is && xid < ie && yid >= js && yid < je && zid >= ks && zid < ke) {
        // get quantities from dev_conserved
        d_gas = dev_conserved[id];
        d_dust = dev_conserved[5*n_cells + id];
        E = dev_conserved[4*n_cells + id];
        // make sure thread hasn't crashed
        if (E < 0.0 || E != E) return;

        vx =  dev_conserved[1*n_cells + id] / d_gas;
        vy =  dev_conserved[2*n_cells + id] / d_gas;
        vz =  dev_conserved[3*n_cells + id] / d_gas;

        #ifdef DE
        ge = dev_conserved[(n_fields-1)*n_cells + id] / d_gas;
        ge = fmax(ge, (double) TINY_NUMBER);
        #endif // DE

        // calculate physical quantities
        p = Calc_Pressure(E, d_gas, vx, vy, vz, gamma);

        double T_init;
        T_init = Calc_Temp(p, n);

        #ifdef DE
        T_init = Calc_Temp_DE(d_gas, ge, gamma, n);
        #endif // DE

        T = T_init;

        // calculate change in dust density
        Dust dustObj(T, n, dt, d_gas, d_dust);
        dustObj.calc_tau_sp();

        dd_dt = dustObj.calc_dd_dt();
        dd = dd_dt * dt;

        // ensure that dust density is not changing too rapidly
        while (d_dust/dd > dd_max) {
            dt_sub = dd_max * d_dust / dd_dt;
            dustObj.d_dust += dt_sub * dd_dt;
            dustObj.dt -= dt_sub;
            dt = dustObj.dt;
            dd_dt = dustObj.calc_dd_dt();
            dd = dt * dd_dt;
        }

        // update dust and gas densities
        dev_conserved[5*n_cells + id] = dustObj.d_dust;
        dev_conserved[id] += dd;
    }

}

void Dust::calc_tau_sp() {
  double a1 = 1; // dust grain size in units of 0.1 micrometers
  double d0 = n/(6*pow(10,-4)); // gas density in units of 10^-27 g/cm^3
  double T_0 = 2*pow(10,6); // K
  double omega = 2.5;
  double A = 0.17*pow(10,9) * YR_IN_S; // 0.17 Gyr in s

  tau_sp = A * (a1/d0) * (pow(T_0/T, omega) + 1); // s
}

double Dust::calc_dd_dt() {
    return -d_dust / (tau_sp/3);
}

void Get_Indices(int n_ghost, int nx, int ny, int nz, int &is, int &ie, int &js, int &je, int &ks, int &ke) {
    is = n_ghost;
    ie = nx - n_ghost;
    if (ny == 1) {
    js = 0;
    je = 1;
    } else {
    js = n_ghost;
    je = ny - n_ghost;
    }
    if (nz == 1) {
    ks = 0;
    ke = 1;
    } else {
    ks = n_ghost;
    ke = nz - n_ghost;
    }
}

void Get_GTID(int &id, int &xid, int &yid, int &zid, int &tid, int nx, int ny, int nz) {
    int blockId = blockIdx + blockIdy * gridDimx;
    id = threadIdx + blockId * blockDimx;
    zid = id / (nx * ny);
    yid = (id - zid * nx * ny) / nx;
    xid = id - zid * nx * ny - yid * nx;
    // add a thread id within the block
    tid = threadIdx;
}

double Calc_Pressure(double E, double d_gas, double vx, double vy, double vz, double gamma) {
    double p;
    p  = (E - 0.5 * d_gas * (vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    p  = fmax(p, (double) TINY_NUMBER);
    return p;
}

double Calc_Temp(double p, double n) {
    double T = p * PRESSURE_UNIT / (n * KB);
    return T;
}

double Calc_Temp_DE(double d_gas, double ge, double gamma, double n) {
    double T =  d_gas * ge * (gamma - 1.0) * PRESSURE_UNIT / (n * KB);
    return T;
}