#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>

using namespace std;


void save_psi(const vector<vector<double>>& psi, const string& filename){
    double NR = 256;
    double NZ = 512;
    ofstream fout(filename);
    fout << fixed << setprecision(8);
    for (int j = 0; j < NZ; ++j) {
        for (int i = 0; i < NR; ++i) {
            fout << psi[i][j];
            if (i < NR - 1) fout << ",";
        }
        fout << "\n";
    }
    fout.close();
    cout << "Saved: " << filename << endl;
}

double psi_analytic(double r, double z){
    double mu0 = 4 * M_PI * 1e-7;
    double B0  = 0.5;
    double r0 = 0.95;
    double k0 = 1.8;
    double q0 = 1.1;
    double term1 = 0.5 * (B0 / (r0*r0 * k0 * q0)) * r*r * z*z;
    double term2 = 0.125 * (B0 * k0*k0 / (r0 * k0 * q0)) * pow((r*r - r0*r0),2);
    return term1 + term2;

}


void enforce_boundaries(std::vector<std::vector<double>>& psi,
    double R_min, double R_max,
    double Z_min, double Z_max,
    double dR, double dZ,
    int NR, int NZ) {
        
    // Enforce R-boundaries (left and right)
    for (int j = 2; j < NZ + 2; ++j) {
        double z = Z_min + (j-2) * dZ;

        // Left boundary (r = R_min)
        psi[2][j] = psi_analytic(R_min, z);
        psi[3][j] = psi_analytic(R_min + dR, z);
        psi[1][j] = 2 * psi[2][j] - psi[3][j];   // Ghost cell i=1
        psi[0][j] = 3 * psi[2][j] - 2 * psi[3][j]; // Ghost cell i=0

        // Right boundary (r = R_max)
        psi[NR+1][j] = psi_analytic(R_max, z);
        psi[NR][j] = psi_analytic(R_max - dR, z);
        psi[NR+2][j] = 2 * psi[NR+1][j] - psi[NR][j];   // Ghost cell i=NR+2
        psi[NR+3][j] = 3 * psi[NR+1][j] - 2 * psi[NR][j]; // Ghost cell i=NR+3
    }

    // Enforce Z-boundaries (bottom and top)
    for (int i = 2; i < NR + 2; ++i) {
        double r = R_min + (i-2) * dR;

        // Bottom boundary (z = Z_min)
        psi[i][2] = psi_analytic(r, Z_min);
        psi[i][3] = psi_analytic(r, Z_min + dZ);
        psi[i][1] = 2 * psi[i][2] - psi[i][3];   // Ghost cell j=1
        psi[i][0] = 3 * psi[i][2] - 2 * psi[i][3]; // Ghost cell j=0

        // Top boundary (z = Z_max)
        psi[i][NZ+1] = psi_analytic(r, Z_max);
        psi[i][NZ] = psi_analytic(r, Z_max - dZ);
        psi[i][NZ+2] = 2 * psi[i][NZ+1] - psi[i][NZ];   // Ghost cell j=NZ+2
        psi[i][NZ+3] = 3 * psi[i][NZ+1] - 2 * psi[i][NZ]; // Ghost cell j=NZ+3
    }
}

void solve_4th_order() {
    int NR = 256;          // Grid points in R
    int NZ = 512;          // Grid points in Z
    double R_min = 0.15;   // Physical R-domain [R_min, R_max]
    double R_max = 1.5;
    double Z_min = -1.35;  // Physical Z-domain [Z_min, Z_max]
    double Z_max = 1.35;
    double dR = (R_max - R_min) / (NR - 1);
    double dZ = (Z_max - Z_min) / (NZ - 1);

    double tolerance = 1e-8;  // Stop iterating when this difference is reached between two successive iterations
    int max_iter = 50000;
    double omega = 1.5;    // Relaxation factor

    // Solov'ev constants
    const double mu0 = 4 * M_PI * 1e-7;
    const double B0  = 0.5;
    const double r0 = 0.95;
    const double k0 = 1.8;
    const double q0 = 1.1;

    const double a = -(B0 * (k0*k0 + 1)) / (mu0 * r0*r0 * k0 * q0); // dp/dψ
    const double b = 0;  // F dF/dψ

    std::vector<std::vector<double>> psi(NR + 4, std::vector<double>(NZ + 4, 0.0));
    std::vector<std::vector<double>> rhs(NR, std::vector<double>(NZ, 0.0));

    // Initialize RHS
    for (int i = 0; i < NR; ++i) {
    double r = R_min + i * dR;
        for (int j = 0; j < NZ; ++j) {
            double z = Z_min + j * dZ;
            rhs[i][j] = -mu0 * r * r * a - b;
            //rhs[i][j] = r*r*B0*(k0*k0 + 1) / (r0*r0*k0*q0);
        }
    }

    // Initialize psi 
    for (int i = 2; i < NR + 4; ++i) {
        for (int j = 2; j < NZ + 4; ++j) {
            psi[i][j] = 0;
        }
    }

    // Main SOR iteration loop
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;
        double sum_diff = 0.0;
        enforce_boundaries(psi, R_min, R_max, Z_min, Z_max, dR, dZ, NR, NZ);
        for (int i = 2; i < NR + 2; ++i) {
            double r = R_min + (i - 2) * dR;
                for (int j = 2; j < NZ + 2; ++j) {
                double psi_old = psi[i][j];

                // 4th-order calculation
                double d2psi_dr2 = (-psi[i+2][j] + 16*psi[i+1][j] - 30*psi[i][j] + 16*psi[i-1][j] - psi[i-2][j]) / (12 * dR * dR);
                double d2psi_dz2 = (-psi[i][j+2] + 16*psi[i][j+1] - 30*psi[i][j] + 16*psi[i][j-1] - psi[i][j-2]) / (12 * dZ * dZ);

                double dpsi_dr;
                dpsi_dr = (-psi[i+2][j] + 8*psi[i+1][j] - 8*psi[i-1][j] + psi[i-2][j]) / (12 * r * dR);
                double test = (rhs[i-2][j-2] - (-psi[i+2][j] + 16*psi[i+1][j]+ 16*psi[i-1][j] - psi[i-2][j])/(12 * dR * dR)\
                            - (-psi[i][j+2] + 16*psi[i][j+1] + 16*psi[i][j-1] - psi[i][j-2]) / (12 * dZ * dZ) + \
                            dpsi_dr) / (-30.0/(12*dR*dR) - 30.0/(12*dZ*dZ)) ;

                double Delta_star_psi = d2psi_dr2 - dpsi_dr / r + d2psi_dz2;
                // Stable SOR update
                double diagonal_weight = -30.0/(12*dR*dR) - 30.0/(12*dZ*dZ);
                double residual = rhs[i-2][j-2] - Delta_star_psi;
                double psi_new = psi_old + omega * residual / diagonal_weight;
                
                double psi_test = (1-omega)*psi_old + omega*test;
                

                psi[i][j] = psi_test;
                double current_diff = std::abs(psi_new - psi_old);
                if (std::isnan(current_diff)) {
                std::cerr << "NaN detected at (" << i << "," << j << ")" << std::endl;
                exit(1);
                }
                max_diff = std::max(max_diff, current_diff);
                sum_diff += current_diff;
            }
        }

        enforce_boundaries(psi, R_min, R_max, Z_min, Z_max, dR, dZ, NR, NZ);

            if (iter % 100 == 0) {
            std::cout << "Iter " << iter << ", max_diff = " << max_diff 
            << ", avg_diff = " << sum_diff/(NR*NZ)
            << ", center = " << psi[NR/2+2][NZ/2+2] << std::endl;
            }
            if (max_diff < tolerance) {
                std::cout << "Converged after " << iter << " iterations." << std::endl;
                break;
            }
    }
    // --- L1 Norm Calculation ---
    double L1_error = 0.0;

    for (int i = 0; i < NR; ++i) {
        double r = R_min + i * dR;
        for (int j = 0; j < NZ; ++j) {
            double z = Z_min + j * dZ;
            double psi_exact = psi_analytic(r, z);
            L1_error += std::abs(psi[i][j] - psi_exact);
        }
    }

    L1_error = L1_error*dR*dZ;
    std::cout << "L1 Norm of error: " << L1_error << std::endl;

    // Extract physical domain for saving
    std::vector<std::vector<double>> psi_physical(NR, std::vector<double>(NZ, 0.0));
    for (int i = 0; i < NR; ++i) {
        for (int j = 0; j < NZ; ++j) {
            psi_physical[i][j] = psi[i+2][j+2];
        }
    }
    save_psi(psi_physical, "gs_solution.csv");
}


int main(){
    auto start = std::chrono::high_resolution_clock::now();
    solve_4th_order();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    return 0;
}

