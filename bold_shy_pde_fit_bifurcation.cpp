#define _USE_MATH_DEFINES

/* cd "/Users/tianxuwang/Documents/Research/bold-shy pde/code" */
/* g++ -std=c++11 -o bold_shy_pde_fit_bifurcation bold_shy_pde_fit_bifurcation.cpp */
/* ./bold_shy_pde_fit_bifurcation */

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>   // for min/max (optional, but we'll use loops)

using std::cin;
using std::cout;
using std::endl;

// ------------------------------------------------------------------
// Memory management templates (unchanged)
// ------------------------------------------------------------------
template <typename T>
T *new_1D(int Nx)
{
    T *ptr = new T[Nx];
    return ptr;
}

template <typename T>
void delete_1D(T *ptr)
{
    delete[] ptr;
}

template <typename T>
T **new_2D(int Nx, int Ny)
{
    T **ptr = new T* [Nx];
    size_t Nxy = (size_t)Nx * Ny;
    ptr[0] = new T[Nxy];

    for(int i = 1; i < Nx; ++i)
        ptr[i] = ptr[i - 1] + Ny;
    return ptr;
}

template <typename T>
void delete_2D(T **ptr)
{
    delete[] ptr[0];
    delete[] ptr;
}

// ------------------------------------------------------------------
// Global parameters (most are constexpr, except xi1 and xi2)
// ------------------------------------------------------------------
constexpr double r = 1.2;
constexpr double K = 1.5;
constexpr double e1 = 0.8;
constexpr double e2 = 0.5;
constexpr double delta1 = 1;
constexpr double delta2 = 1;
constexpr double delta3 = 1;

constexpr double a1 = -0;
constexpr double b1 = 0;
constexpr double c1 = 6;
constexpr double a2 = -0.;
constexpr double b2 = 0.;
constexpr double c2 = 0.07;
// xi1 and xi2 will be varied – make them non‑const global variables
double xi1;
double xi2;
constexpr double eta = 1;
constexpr double tau1 = 1.235;
constexpr double tau2 = 33.33;

constexpr double d11 = 0.001;
constexpr double d20 = 0.25;
constexpr double d30 = 0.003;
constexpr double d2 = 0.0159;
constexpr double d31 = 0.04;
constexpr double beta1 = 0.00;
constexpr double beta2 = 0.00;
constexpr double T0 = 0;
constexpr double kappa = 0.1;

constexpr double A = 100;
constexpr double B = 15;
constexpr int    Nx = 1000;
constexpr int    Nt = 300000;
constexpr int    sample_itv = 100;
constexpr double T_final = 2500;
constexpr double dt = T_final / (double)Nt;
constexpr double dx = 2.0 * A / (double)Nx;

// ------------------------------------------------------------------
// Functions that depend on temperature (which is zero here, but kept general)
// ------------------------------------------------------------------
static double T(double x) { return T0 * exp(-x * x / (B * B)); }

static double mu11(double T) { return fmax(a1 * T * T + b1 * T + c1, c1); }
static double mu21(double T) { return fmax(a2 * T * T + b2 * T + c2, c2); }

// These now use the global variables xi1, xi2
static double mu12(double T) { return xi1 * mu11(T); }
static double mu22(double T) { return xi2 * mu21(T); }

static double alpha1(double T) { return eta * mu11(T); }
static double alpha2(double T) { return eta * mu12(T); }
static double gamma1(double T) { return eta * mu21(T); }
static double gamma2(double T) { return eta * mu22(T); }

static double d1(double T)   { return d11 * T; }
static double d21(double T)  { return d20 + d2 * T; }
static double d22(double T)  { return d20 + d2 * T; }
static double d3(double T)   { return d30 + d31 * T; }

static double h(double z) { return z / (1 + kappa * z); }

// Initial conditions (small perturbation)
static double X0 (double x) { return 0.5 + 0.005 * cos(x / 100 * 7); }
static double Y10(double x) { return 0.5 - 0.005 * cos(x / 100 * 5); }
static double Y20(double x) { return 0.5 - 0.005 * cos(x / 100 * 5); }
static double Z0 (double x) { return 10 + 0.005 * cos(x / 100 * 9); }

// ------------------------------------------------------------------
// Interpolation and source terms (unchanged)
// ------------------------------------------------------------------
void interpolate_half(double **uh, double **uh_half)
{
    for(int i = 0; i < 4; ++i)
    {
        for(int j = 0; j < Nx + 2; ++j)
        {
            if(j == 0 || j == Nx + 1)
                uh_half[i][j] = (uh[i][j] + uh[i][j + 1]) / 2.0;
            else
                uh_half[i][j] = (9.0 * uh[i][j] + 9.0 * uh[i][j + 1] - uh[i][j - 1] - uh[i][j + 2]) / 16.0;
        }
    }
}

void setSourceVector(double **uh, double *mesh, double **source)
{
    for(int j = 0; j < Nx + 3; ++j)
    {
        if(j == 0 || j == Nx + 2)
            for(int i = 0; i < 4; ++i)
                source[i][j] = 0;
        else
        {
            double Tval = T(mesh[j]);
            double fxy1 = mu11(Tval) * uh[0][j] * uh[1][j] / (1.0 + mu11(Tval) * tau1 * uh[0][j]);
            double fxy2 = mu12(Tval) * uh[0][j] * uh[2][j] / (1.0 + mu12(Tval) * tau1 * uh[0][j]);
            double gy1z = mu21(Tval) * uh[1][j] * uh[3][j] / (1.0 + mu21(Tval) * tau2 * uh[1][j]);
            double gy2z = mu22(Tval) * uh[2][j] * uh[3][j] / (1.0 + mu22(Tval) * tau2 * uh[2][j]);

            source[0][j] = r * uh[0][j] * (1 - uh[0][j] / K) - fxy1 - fxy2 - d1(Tval) * uh[0][j];
            source[1][j] = e1 * fxy1 - gy1z - d21(Tval) * uh[1][j] - (beta1 * Tval) * uh[1][j] + (beta2 * Tval) * uh[2][j];
            source[2][j] = e1 * fxy2 - gy2z - d22(Tval) * uh[2][j] + (beta1 * Tval) * uh[1][j] - (beta2 * Tval) * uh[2][j];
            source[3][j] = e2 * gy1z + e2 * gy2z - d3(Tval) * uh[3][j];
        }
    }
}

// ------------------------------------------------------------------
// Main function: two‑parameter scan over xi1 and xi2
// ------------------------------------------------------------------
int main()
{
    // ------------------------------------------------------------------
    // Parameter ranges for the bifurcation diagram
    // ------------------------------------------------------------------
    const double xi1_start = 0.1;
    const double xi1_end   = 1.;
    const double xi1_step  = 0.05;

    const double xi2_start = 0.1;
    const double xi2_end   = 1.;
    const double xi2_step  = 0.05;

    const int num_xi1 = static_cast<int>((xi1_end - xi1_start) / xi1_step) + 1;
    const int num_xi2 = static_cast<int>((xi2_end - xi2_start) / xi2_step) + 1;

    // Spatial point where we monitor Y1 and Y2 (centre of domain)
    const int spatial_index = Nx / 2;          // 0‑based index among interior points
    const int j_center = spatial_index + 1;    // corresponding index in uh (ghost cells at 0 and Nx+2)

    // Fraction of the simulation used for the "last short period"
    const double last_fraction = 0.25;           // use last 30%

    // Threshold for survival (average maximum must be above this value)
    const double survival_threshold = 5e-3;

    // ------------------------------------------------------------------
    // Pre‑compute mesh (independent of xi1, xi2)
    // ------------------------------------------------------------------
    double *mesh = new_1D<double>(Nx + 3);
    double *half_mesh = new_1D<double>(Nx + 2);

    for(int j = -1; j <= Nx + 1; ++j)
        mesh[j + 1] = j * dx - A;
    for(int j = 0; j <= Nx + 1; ++j)
        half_mesh[j] = (j - 0.5) * dx - A;

    // ------------------------------------------------------------------
    // Open output file (write header)
    // ------------------------------------------------------------------
    std::ofstream fout("bifurcation_2D.txt", std::ios::trunc);
    fout << "# xi1 xi2 survive_Y1 survive_Y2 amplitude_Y1 amplitude_Y2\n";
    fout.close();

    // ------------------------------------------------------------------
    // Outer loop over xi1
    // ------------------------------------------------------------------
    for(int i1 = 0; i1 < num_xi1; ++i1)
    {
        xi1 = xi1_start + i1 * xi1_step;
        cout << "xi1 = " << xi1 << " : " << std::flush;

        // Inner loop over xi2
        for(int i2 = 0; i2 < num_xi2; ++i2)
        {
            xi2 = xi2_start + i2 * xi2_step;
            cout << xi2 << " " << std::flush;

            clock_t run_begin = clock();

            // ------------------------------------------------------------------
            // Allocate per‑run arrays (depend on xi1, xi2 via alpha2, gamma2 etc.)
            // ------------------------------------------------------------------
            double **uh      = new_2D<double>(4, Nx + 3);
            double **uh_half = new_2D<double>(4, Nx + 2);
            double **source  = new_2D<double>(4, Nx + 3);
            double **temp    = new_2D<double>(4, Nx + 3);

            // Coefficients on half‑mesh (depend on xi1, xi2 through mu12, mu22)
            double *alpha1_half = new_1D<double>(Nx + 2);
            double *alpha2_half = new_1D<double>(Nx + 2);
            double *gamma1_half = new_1D<double>(Nx + 2);
            double *gamma2_half = new_1D<double>(Nx + 2);

            for(int j = 0; j <= Nx + 1; ++j)
            {
                double T_val = T(half_mesh[j]);
                alpha1_half[j] = alpha1(T_val);
                alpha2_half[j] = alpha2(T_val);
                gamma1_half[j] = gamma1(T_val);
                gamma2_half[j] = gamma2(T_val);
            }

            // ------------------------------------------------------------------
            // Set initial conditions
            // ------------------------------------------------------------------
            for(int j = 0; j < Nx + 3; ++j)
            {
                uh[0][j] = X0 (mesh[j]);
                uh[1][j] = Y10(mesh[j]);
                uh[2][j] = Y20(mesh[j]);
                uh[3][j] = Z0 (mesh[j]);
            }

            // ------------------------------------------------------------------
            // Vectors to store time series of Y1 and Y2 at the centre
            // ------------------------------------------------------------------
            std::vector<double> series_Y1, series_Y2;
            // Record initial state
            series_Y1.push_back(uh[1][j_center]);
            series_Y2.push_back(uh[2][j_center]);

            // ------------------------------------------------------------------
            // Time integration
            // ------------------------------------------------------------------
            int iter = 0;
            while(iter < Nt)
            {
                interpolate_half(uh, uh_half);
                setSourceVector(uh, mesh, source);

                for(int j = 1; j <= Nx + 1; ++j)
                {
                    temp[0][j] = dt * source[0][j] + uh[0][j] + dt / dx / dx * delta1 * (uh[0][j - 1] + uh[0][j + 1] - 2.0 * uh[0][j]);
                    temp[1][j] = dt * source[1][j] + uh[1][j] + dt / dx / dx *
                        ( delta2 * (uh[1][j - 1] + uh[1][j + 1] - 2.0 * uh[1][j])
                        + alpha1_half[j - 1] * h(uh_half[1][j - 1]) * (uh[0][j] - uh[0][j - 1])
                        + alpha1_half[j]     * h(uh_half[1][j])     * (uh[0][j] - uh[0][j + 1]) );
                    temp[2][j] = dt * source[2][j] + uh[2][j] + dt / dx / dx *
                        ( delta2 * (uh[2][j - 1] + uh[2][j + 1] - 2.0 * uh[2][j])
                        + alpha2_half[j - 1] * h(uh_half[2][j - 1]) * (uh[0][j] - uh[0][j - 1])
                        + alpha2_half[j]     * h(uh_half[2][j])     * (uh[0][j] - uh[0][j + 1]) );
                    temp[3][j] = dt * source[3][j] + uh[3][j] + dt / dx / dx *
                        ( delta3 * (uh[3][j - 1] + uh[3][j + 1] - 2.0 * uh[3][j])
                        + gamma1_half[j - 1]  * h(uh_half[3][j - 1]) * (uh[1][j] - uh[1][j - 1])
                        + gamma1_half[j]      * h(uh_half[3][j])     * (uh[1][j] - uh[1][j + 1])
                        + gamma2_half[j - 1]  * h(uh_half[3][j - 1]) * (uh[2][j] - uh[2][j - 1])
                        + gamma2_half[j]      * h(uh_half[3][j])     * (uh[2][j] - uh[2][j + 1]) );
                }

                for(int i = 0; i < 4; ++i)
                    for(int j = 1; j <= Nx + 1; ++j)
                        uh[i][j] = fmax(0., temp[i][j]);
                for(int i = 0; i < 4; ++i)
                {
                    uh[i][0] = uh[i][2];
                    uh[i][Nx + 2] = uh[i][Nx];
                }

                ++iter;
                if(iter % sample_itv == 0)
                {
                    // Record current values at the centre
                    series_Y1.push_back(uh[1][j_center]);
                    series_Y2.push_back(uh[2][j_center]);
                }
            }

            // ------------------------------------------------------------------
            // Compute statistics from the last part of the simulation
            // ------------------------------------------------------------------
            int total_rec = series_Y1.size();   // same for Y2
            int last_start = total_rec - static_cast<int>(total_rec * last_fraction);
            if (last_start < 0) last_start = 0;
            if (last_start > total_rec - 3) last_start = total_rec - 3;   // need at least three points

            // --- Average of maxima (survival criterion) ---
            auto average_max = [&](const std::vector<double>& series) -> double {
                double sum = 0.0;
                int count = 0;
                for (size_t i = last_start + 1; i + 1 < series.size(); ++i) {
                    if (series[i-1] < series[i] && series[i] > series[i+1]) {
                        sum += series[i];
                        ++count;
                    }
                }
                return (count > 0) ? (sum / count) : 0.0;
            };

            double avg_max_Y1 = average_max(series_Y1);
            double avg_max_Y2 = average_max(series_Y2);

            int survive_Y1 = (avg_max_Y1 > survival_threshold) ? 1 : 0;
            int survive_Y2 = (avg_max_Y2 > survival_threshold) ? 1 : 0;

            // --- Minimum and maximum (to compute amplitude) ---
            double min_Y1 = series_Y1[last_start];
            double max_Y1 = series_Y1[last_start];
            for (int i = last_start + 1; i < total_rec; ++i) {
                if (series_Y1[i] < min_Y1) min_Y1 = series_Y1[i];
                if (series_Y1[i] > max_Y1) max_Y1 = series_Y1[i];
            }
            double amplitude_Y1 = max_Y1 - min_Y1;

            double min_Y2 = series_Y2[last_start];
            double max_Y2 = series_Y2[last_start];
            for (int i = last_start + 1; i < total_rec; ++i) {
                if (series_Y2[i] < min_Y2) min_Y2 = series_Y2[i];
                if (series_Y2[i] > max_Y2) max_Y2 = series_Y2[i];
            }
            double amplitude_Y2 = max_Y2 - min_Y2;

            // ------------------------------------------------------------------
            // Append result to output file
            // ------------------------------------------------------------------
            fout.open("bifurcation_2D.txt", std::ios::app);
            fout << xi1 << " " << xi2 << " "
                 << survive_Y1 << " " << survive_Y2 << " "
                 << amplitude_Y1 << " " << amplitude_Y2 << "\n";
            fout.close();

            // ------------------------------------------------------------------
            // Clean up per‑run arrays
            // ------------------------------------------------------------------
            delete_2D(uh);
            delete_2D(uh_half);
            delete_2D(source);
            delete_2D(temp);
            delete_1D(alpha1_half);
            delete_1D(alpha2_half);
            delete_1D(gamma1_half);
            delete_1D(gamma2_half);
        }
        cout << endl;
    }

    // ------------------------------------------------------------------
    // Final clean‑up
    // ------------------------------------------------------------------
    delete_1D(mesh);
    delete_1D(half_mesh);

    cout << "All done." << endl;
    return 0;
}
