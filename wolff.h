#ifndef WOLFF_H
#define WOLFF_H

#include <iostream>
#include <complex>
#include <cmath>
#include <string>
#include <filesystem>
#include <fstream>

using std::complex;
using std::conj;
using std::string;

class WOLFF {
    private:
        int q;
        int L;
        int N;
        int *spins;

        int *m_counter;
        complex<double> *m_weights;

        double last_T;
        complex<double> *corr_function;
        complex<double> m;
        double m1, m2, m4;

        bool write_header = true;

        enum dirs {RIGHT, LEFT};

        int indx(int x) {return x;}
        int xpos(int i) {return i % L;}

        int nbr(int i, int dir) {
            int x = xpos(i);

            switch (dir) {
                case RIGHT: return indx((x + 1) % L);
                case LEFT: return indx((x - 1 + L) % L);
            }
            return -1;
        }

        void flip_and_build_from(int s, double T) {
            int oldstate = spins[s], newstate = (spins[s] + 1) % q;
            
            spins[s] = newstate;
            m_counter[oldstate]--;
            m_counter[newstate]++;

            for (int dir = 0; dir < 2; dir++) {
                int j = nbr(s, dir);
                if (spins[j] == oldstate && rand() / (RAND_MAX + 1.) < p_connect(T)) {
                    flip_and_build_from(j, T);
                }
            }
        }

        double p_connect(double T) {
            return 1 - exp(- 1 / T);
        }

    public:

    WOLFF(int q_, int L_, int N_) {
        int i, s;
        srand((unsigned) time(0));

        q = q_;
        L = L_;
        N = N_;

        spins = new int[N];

        m_counter = new int[q];
        m_weights = new complex<double>[q];

        corr_function = new complex<double>[N + 1];

        for (s = 0; s < q; s++) {
            m_weights[s] = complex<double>(cos(2 * M_PI * s / q), sin(2 * M_PI * s / q));
            m_counter[s] = 0;
        }
        m_counter[0] = N;
        
        for (i = 0; i < N; i++) {
            spins[i] = 0;
        }
    }

    ~WOLFF() {
        delete[] spins;
        delete[] m_counter;
        delete[] m_weights;
        delete[] corr_function;
    }

    void simulate(double T, long therm_steps, long mc_cycles, int n_bins = 10, int n_clusters = 1) {
        last_T = T;

        for (int t = 0; t < therm_steps; t++) {
            for (int c = 0; c < n_clusters; c++) {
                flip_and_build_from(rand() % N, T);
            }
        }
            
        complex<double> m0;
        complex<double> *mr = new complex<double>[N];
        complex<double> *m0r = new complex<double>[N];
        for (int r = 0; r < N; r++) {
            corr_function[r] = 0.0;
        }
        m = 0.0;
        m1 = 0.0;
        m2 = 0.0;
        m4 = 0.0;            

        for (int n = 0; n < n_bins; n++) {
            m0 = 0.0;
            for (int i = 0; i < N; i++) {
                mr[i] = 0.0;
                m0r[i] = 0.0;
            }

            for (int t = 0; t < mc_cycles; t++) {
                for (int c = 0; c < n_clusters; c++) {
                    flip_and_build_from(rand() % N, T);
                }
                
                complex<double> tm(0.0, 0.0);
                for (int s = 0; s < q; s++) {
                    tm += m_weights[s] * (double) m_counter[s];
                }
                tm /= N;
                double tm1 = abs(tm);
                double tm2 = tm1 * tm1;
                m += tm; m1 += tm1; m2 += tm2; m4 += tm2 * tm2;

                m0 += conj(m_weights[spins[0]]);
                for (int r = 0; r < N; r++) {
                    mr[r] += m_weights[spins[r]];
                    m0r[r] += conj(m_weights[spins[0]]) * m_weights[spins[r]];
                }
            }

            m0 /= mc_cycles;
            for (int r = 0; r < N; r++) {
                mr[r] /= mc_cycles;
                m0r[r] /= mc_cycles;

                corr_function[r] += m0r[r] - m0 * mr[r];
            }
        }

        for (int r = 0; r < N; r++) {
            corr_function[r] /= n_bins;
        }

        m /= (mc_cycles * n_bins); 
        m1 /= (mc_cycles * n_bins); 
        m2 /= (mc_cycles * n_bins); 
        m4 /= (mc_cycles * n_bins);

        delete[] mr;
        delete[] m0r;
    }

    void print_results() {
        printf("T = %f \n", last_T);
        printf("<m> = %.5f | <|m|> = %.5f | <|m2|> = %.5f | <|m4|> = %.5f \n", m.real(), m1, m2, m4);
        for (int r = 0; r < N + 1; r++) {
            printf("C[%i] = %.5f \n", r, corr_function[r % N].real());
        }
        printf("\n");
    }

    void write_results(string path, string name) {
        std::filesystem::create_directories(path); 

        std::ofstream file(path + name, std::ios_base::app);
        if (file.is_open()) {
            if (write_header) {
                file << "N,L,q\n";
                file << N << "," << L << "," << q << "\n";
                write_header = false;
            }
            file << last_T << "\n";
            file << m.real() << "," << m1 << "," << m2 << "," << m4 << "\n";
            for (int r = 0; r < N + 1; r++) {
                file << r << "," << corr_function[r % N].real() << "\n";
            }
            file.close();
        } else {
            printf(" -- Error: can not open save file, please check you directory -- \n");
        }
    }

};

#endif // WOLFF_H

