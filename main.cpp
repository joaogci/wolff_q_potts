#include <iostream>

#include "wolff.h"

int main(int argc, char **argv) {
    int q = 3;
    int L = 32;

    WOLFF wolff(q, L, L);

    double T_vals[] = {0, 0.25, 0.5, 1.0, 2.0, 4.0};
    for (int t_idx = 0; t_idx < 6; t_idx++) {
        wolff.simulate(T_vals[t_idx], 1e5, 1e6, 1);
        wolff.print_results();
        wolff.write_results("../", "1D_L32_q3_potts.csv");
    }

    return 0;
}
