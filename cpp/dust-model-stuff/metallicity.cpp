#include <iostream>
#include<math.h>
#include<vector>


int main() {
    float MP = 1.6726E-24;  // proton mass in grams
    float KB = 1.380658e-16; // Boltzmann constant
    float LENGTH_UNIT = 3.08567758e21; // kpc in cm
    float TINY_NUMBER = 1.0e-20;
    float MASS_UNIT = 1.98855e33; // solar mass in grams
    float TIME_UNIT = 3.15569e10; // kyr in s

    float n = 1e-2;
    float d_gas = MP * n * pow(LENGTH_UNIT, 3) / MASS_UNIT;
    float d_dust = d_gas / 100;

    /* mass fractions for metals assuming solar metallicity
        O, C, Ni, Si, Mg, Ne, Fe, S */
    float metallicity = 0;
    float d_metal;
    std::vector<float> metals = {0.97, 0.40, 0.096, 0.099, 0.079, 0.058, 0.14, 0.040};

    for (int i = 0; i < metals.size(); i++) {
        metallicity += metals[i];
    }

    //for (std::metals::iterator i = vector.begin(); i != vector.end()); ++i) {
    //    metallicity += i;
    //}

    d_metal = metallicity * d_gas;

    printf("%f\n", metallicity); 
    printf("%f\n", d_metal);
}