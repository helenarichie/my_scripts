import numpy as np

class DustIntegrator:
    MP = 1.6726E-24;  # proton mass in g
    YR_IN_S = 3.154e7 # one year in s
    # use this ^ conversion just because it's easier to think about time in terms of years
    
    # solar abundances (percentage of total mass)
    # O, C, N, Si, Mg, Ne, Fe, S
    metals = np.array([0.0097, 0.0040, 0.00096, 0.00099, 0.00076, 0.00058, 0.0014, 0.00040]);
    metallicity = np.sum(metals)
    
    def __init__(self, T, n, dt, tmax, int_type):
        self.T = T # K
        self.n = n
        self.dt = dt # passed in in s
        self.tmax = tmax # passed in in s
        self._d0_gas, self._d0_metal, self._d0_dust = self.calc_init_density() # g/cm^3
        self.d_gas, self.d_metal, self.d_dust = self.calc_init_density() # g/cm^3
        self.tau_g = self.calc_tau_g() # s
        self.tau_sp = self.calc_tau_sp() # s
        self.int_type = int_type # tot, acc, or sp

    def calc_tau_g(self):
        tau_g_ref = 0.2e9 * self.YR_IN_S # 0.2 Gyr in s
        d_ref = self.MP # 1 H atom per cubic centimer
        T_ref = 20.0; # 20 K
        tau_g = tau_g_ref * (d_ref/(self._d0_gas)) * (T_ref/self.T) ** (1/2); # s
    
        return tau_g
    
    def calc_tau_sp(self):
        a1 = 1; # dust grain size in units of 0.1 micrometer
        d0 = self.n/(6e-4); # gas density in units of 10^-27 g/cm^3
        T_0 = 2e6; # K
        omega = 2.5;
        A = 0.17e9 * self.YR_IN_S; # 0.17 Gyr in s
        
        tau_sp = A * (a1/d0) * ((T_0/self.T) ** omega + 1); # s
        
        return tau_sp
    
    def calc_init_density(self):
        d0_gas = self.MP * self.n # g/cm^3
        d0_metal = self.metallicity * d0_gas
        # assume 1% dust to gas fraction
        d0_dust = d0_gas / 100
        
        return d0_gas, d0_metal, d0_dust

    def calc_dd_dt(self):
        # calculate the rate of change of dust density
        # g/cm^3/s
        dd_dt = None
        if self.int_type == "tot":
            dd_dt = (1 - self.d_dust/self.d_metal)*(self.d_dust/self.tau_g) - self.d_dust/(self.tau_sp/3)
        if self.int_type == "acc":
            dd_dt = (1 - self.d_dust/self.d_metal)*(self.d_dust/self.tau_g)
        if self.int_type == "sp":
            dd_dt = - self.d_dust/(self.tau_sp/3)
        
        return dd_dt

    def calc_d_dust(self):
        # calculate the new total dust density
        dd_dt = self.calc_dd_dt()
        d_dust_i = self.d_dust + self.dt * dd_dt
        
        return d_dust_i
    
    def update_densities(self):
        # update the integrator's density attributes
        dd_dt = self.calc_dd_dt()
        self.d_dust += self.dt * dd_dt
        if self.int_type == "acc":
            self.d_metal -= self.dt * dd_dt
    
    def get_info(self):
        print("Gas temperature: {:.1E} K".format(self.T))
        print("Gas initial mass density: {:.5E} g/cm^3".format(self._d0_gas))
        print("Gas initial number density: {} cm^-3".format(self.n))
        print("Growth timescale: {:.5E} yr".format(self.tau_g/self.YR_IN_S))
        print("Destruction timescale: {:.5E} yr".format(self.tau_sp/self.YR_IN_S))
        print("Time-step: {} yr".format(self.dt/self.YR_IN_S))
        print("------------------------------------------")

def method():
    pass

def evolve_solutions(h, n, T_arr, tmax, int_type):
    YEAR_IN_SECONDS = 3.154e7
    h *= YEAR_IN_SECONDS
    tmax *= YEAR_IN_SECONDS
    
    t_arr = np.arange(0, tmax, h)
    
    d_dust = np.zeros(shape=(len(T_arr), len(t_arr)))
    d_metal = np.zeros(shape=(len(T_arr), len(t_arr)))

    integrators = []
    
    tau_gs = []
    tau_sps = []
    
    def time_refine(integrator):
        """Refine time-step if dd is changing too rapidly."""
        d_dust_i = integrator.calc_d_dust()
        dd_dt = None
        if self.int_type == "tot":
            dd_dt = self.calc_dd_dt_tot()
        if self.int_type == "acc":
            dd_dt = self.calc_dd_dt_acc()
        if self.int_type == "sp":
            dd_dt = self.calc_dd_dt_sp()
        dt = integrator.dt
        dd = dd_dt * dt
        while dd/d_dust_i > 0.01:
            dt_sub = 0.01 * d_dust_i / dd_dt

            integrator.d_dust += dt_sub * dd_dt
            if integrator.int_type == "acc":
                integrator.d_metal -= dt_sub * dd_dt
            integrator.dt -= dt_sub
            dt = integrator.dt
            dd_dt = integrator.calc_dd_dt()
            dd = dt * dd_dt

    for j, T in enumerate(T_arr):
        integrator = DustIntegrator(T, n, h, tmax, int_type)
        integrator.get_info()
        # set initial densities
        d_dust[j][0] = integrator._d0_dust
        if integrator.int_type == "acc":
            d_metal[j][0] = integrator._d0_metal
        tau_gs.append(integrator.tau_g) # s
        tau_sps.append(integrator.tau_sp) # s
        integrators.append(integrator)

        for i, t_i in enumerate(t_arr):
            if (i+1) < len(t_arr):
                for integrator in integrators:
                    # calculate dust density for this time-step
                    d_dust_i = integrator.calc_d_dust()
                    # calculate dd
                    dd_dt = integrator.calc_dd_dt()
                    dd = integrator.dt * dd_dt

                    # if rate of dust growth is changing too rapidly
                    if dd/d_dust_i > 0.01:
                        # shorten time-step
                        time_refine(integrator)
                        integrator.dt = h
                        d_dust[j][i+1] = integrator.d_dust
                        if integrator.int_type == "acc":
                            d_metal[j][i+1] = integrator.d_metal
                        continue

                    # this is where the density attribute gets updated
                    integrator.update_densities()
                    d_dust[j][i+1] = integrator.d_dust
                    if integrator.int_type == "acc":
                        d_metal[j][i+1] = integrator.d_metal

    return d_dust, d_metal, tau_gs, tau_sps, integrators