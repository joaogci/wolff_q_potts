#include <iostream>
#include <string>

#include "wolff.h"

using std::string;
using std::to_string;

int main(int argc, char **argv) {
    int q = 3;
    int L = 32;
    int dim = 2;
    int N = pow(L, dim);
    string file_name = to_string(dim) + "D_L" + to_string(L) + "_q" + to_string(q) + "_potts.csv";

    WOLFF wolff(q, L, N, dim);

    int size_T = 50;
    double T_s = 0.0;
    double T_e = 2.0;
    double T_vals[size_T];

    for (int t_idx = 0; t_idx < size_T; t_idx++) {
        T_vals[t_idx] = (T_e - T_s) * (t_idx + 1) / size_T;

        wolff.simulate(T_vals[t_idx], 1e5, 1e6, 10);
        wolff.write_results("../", file_name);
        printf("t_idx = %i \n", t_idx);
    }

    return 0;
}
