#include <iostream>

#include "wolff.h"

int main(int argc, char **argv) {
    int q = 3;
    int L = 64;
    int dim = 2;

    WOLFF wolff(q, L, L * L, dim);

    double T_vals[] = {0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};
    for (int t_idx = 0; t_idx < 9; t_idx++) {
        wolff.simulate(T_vals[t_idx], 1e5, 1e6, 1);
        // wolff.print_results();
        wolff.write_results("../", "2D_L64_q3_potts_pt.csv");
        printf("t_idx = %i \n", t_idx);
    }

    return 0;
}
