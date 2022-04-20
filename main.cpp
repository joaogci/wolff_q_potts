#include <iostream>

#include "wolff.h"

int main(int argc, char **argv) {
    int q = 3;
    int L = 16;
    int dim = 2;
    int N = pow(L, dim);

    WOLFF wolff(q, L, N, dim);

    int size_T = 30;
    double T_s = 0;
    double T_e = 3;
    double T_vals[size_T];

    for (int t_idx = 0; t_idx < size_T; t_idx++) {
        T_vals[t_idx] = (T_e - T_s) * t_idx / size_T;

        wolff.simulate(T_vals[t_idx], 1e5, 1e6, 1);
        wolff.write_results("../", "2D_L16_q3_potts.csv");
        printf("t_idx = %i \n", t_idx);
    }

    return 0;
}
