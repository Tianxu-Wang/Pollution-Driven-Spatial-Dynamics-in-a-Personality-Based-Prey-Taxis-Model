#define _USE_MATH_DEFINES

/* cd "/Users/tianxuwang/Documents/Research/bold-shy pde/code" */
/* g++ -std=c++11 -o bold_shy_pde_fit bold_shy_pde_fit.cpp */
/* ./bold_shy_pde_fit */

#include <cmath>
#include <fstream>
#include <iostream>

using std::cin;
using std::cout;
using std::endl;

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

template <typename T>
T ***new_3D(int Nx, int Ny, int Nz)
{
	size_t Nxy  = (size_t)Nx * Ny;
	size_t Nyz  = (size_t)Ny * Nz;
	size_t Nxyz = (size_t)Nx * Ny * Nz;
	
	T ***ptr  = new T** [Nx];
	ptr[0]    = new T* [Nxy];
	ptr[0][0] = new T [Nxyz];

	for(int i = 0; i < Nx; ++i)
	{
		if(i != 0)
		{
			ptr[i] = ptr[i - 1] + Ny;
			ptr[i][0] = ptr[i - 1][0] + Nyz;
		}
		for(int j = 1; j < Ny; ++j)
			ptr[i][j] = ptr[i][j - 1] + Nz;
	}
	return ptr;
}

template <typename T>
void delete_3D(T ***ptr)
{
	delete[] ptr[0][0];
	delete[] ptr[0];
	delete[] ptr;
}

constexpr double r = 1.2;
constexpr double K = 1.5;
constexpr double e1 = 0.8;
constexpr double e2 = 0.5;
constexpr double delta1 = 1;
constexpr double delta2 = 1;
constexpr double delta3 = 1;

constexpr double a1 = -0.0034; /*-0.0034;*/
constexpr double b1 = 0.085; /*0.085;*/
constexpr double c1 = 4.5384; /*4.5384; 3.24*/
constexpr double a2 = -0.0004; /*-0.0004; */
constexpr double b2 = 0.0012; /* 0.0012; */
constexpr double c2 = 0.0516; /*0.0516;*/
constexpr double xi1 = 0.6;
constexpr double xi2 = 0.6;
constexpr double eta = 1;
constexpr double tau1 = 1.235;
constexpr double tau2 = 33.33;

constexpr double d11 = 0.00;
//constexpr double d20 = 0.25;
//constexpr double d30 = 0.003;
constexpr double d2 = 0.0159;
constexpr double d31 = 0.00004;
constexpr double beta1 = 0.0; /*0.5;*/
constexpr double beta2 = 0.0000;
constexpr double T0 = 20;
constexpr double kappa = 0.1;

constexpr double A = 100;
constexpr double B = 15;
constexpr int    Nx = 1000;
constexpr int    Nt = 1000000;
// This is the actual number of grid points in time which is used in the solver.
// should be near or at least Nx^2
constexpr int    sample_itv = 100; 
// We record the solution for each *sample_itv* steps. 
// Nt in Matlab is the same as *Nt / sample_itv* in this program
constexpr double T_final = 2000;
constexpr double dt = T_final / (double)Nt;
constexpr double dx = 2.0 * A / (double)Nx;

static double T(double x) { return T0 * exp(-x * x / (B * B)); }

static double mu11(double T) { return fmax(a1 * T * T + b1 * T + c1, c1); }

static double mu21(double T) { return fmax(a2 * T * T + b2 * T + c2, c2); }

static double mu12(double T) { return xi1 * mu11(T); }

static double mu22(double T) { return xi2 * mu21(T); }

static double alpha1(double T) { return eta * mu11(T); }

static double alpha2(double T) { return eta * mu12(T); }

static double gamma1(double T) { return eta * mu21(T); }

static double gamma2(double T) { return eta * mu22(T); }

static double d1(double T) { return d11 * T; }

static double d21(double T) { return 0.013983 * (exp(0.160884 * T) - 1.0) + 0.25; }

static double d22(double T) { return xi1 * 0.013983 * (exp(0.160884 * T) - 1.0) + 0.25; }

static double d3(double T)  { return 0.010872 * (exp(0.071760 * T) - 1.0) + 0.003; }

static double h(double z) { return z / (1 + kappa * z); }

static double X0 (double x) { return 0.5 + 0.005 * cos(x / 100 * 7); }
static double Y10(double x) { return 0.5 - 0.005 * cos(x / 100 * 5); }
static double Y20(double x) { return 0.5 - 0.005 * cos(x / 100 * 5); }
static double Z0 (double x) { return 0.5 + 0.005 * cos(x / 100 * 9); }

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
			source[1][j] = e1 * fxy1 - gy1z - d21(Tval) * uh[1][j] - (beta1 * Tval + 0.00) * uh[1][j] + (beta2 * Tval + 0.00) * uh[2][j];
			source[2][j] = e1 * fxy2 - gy2z - d22(Tval) * uh[2][j] + (beta1 * Tval + 0.00) * uh[1][j] - (beta2 * Tval + 0.00) * uh[2][j];
			source[3][j] = e2 * gy1z + e2 * gy2z - d3(Tval) * uh[3][j];
		}
	}
}

int main()
{
	clock_t clock_begin = clock();
	cout << "Solving PDE..." << std::flush;

	double *mesh        = new_1D<double> (Nx + 3);
	double *half_mesh   = new_1D<double> (Nx + 2);
	double *alpha1_half = new_1D<double> (Nx + 2);
	double *alpha2_half = new_1D<double> (Nx + 2);
	double *gamma1_half = new_1D<double> (Nx + 2);
	double *gamma2_half = new_1D<double> (Nx + 2);

	for(int j = -1; j <= Nx + 1; ++j)
		mesh[j + 1] = j * dx - A;
	for(int j = 0; j <= Nx + 1; ++j)
	{
		half_mesh[j] = (j - 0.5) * dx - A;
		double T_val = T(half_mesh[j]);
		alpha1_half[j] = alpha1(T_val);
		alpha2_half[j] = alpha2(T_val);
		gamma1_half[j] = gamma1(T_val);
        gamma2_half[j] = gamma2(T_val);
	}

	double **uh      = new_2D<double> (4, Nx + 3);
	double **uh_half = new_2D<double> (4, Nx + 2);
	double ***uh_rec = new_3D<double> (4, Nt / sample_itv + 1, Nx + 1);

	for(int j = 0; j < Nx + 3; ++j)
	{
		uh[0][j] = X0 (mesh[j]);
		uh[1][j] = Y10(mesh[j]);
		uh[2][j] = Y20(mesh[j]);
		uh[3][j] = Z0 (mesh[j]);
	}

	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j <= Nx; ++j)
			uh_rec[i][0][j] = uh[i][j + 1];
	}

	int iter = 0;
	double **source = new_2D<double> (4, Nx + 3);
	double **temp     = new_2D<double> (4, Nx + 3);

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
			for(int i = 0; i < 4; ++i)
				for(int j = 0; j <= Nx; ++j)
					uh_rec[i][iter / sample_itv][j] = uh[i][j + 1];
		}
	}
	cout << "  Done." << endl;
	cout << "time elapsed: " << (double)(clock() - clock_begin) / CLOCKS_PER_SEC << "s" << endl;

	cout << "Saving data to file..." << std::flush;
	std::ofstream fout[4];
	fout[0].open("./X.txt",  std::fstream::trunc);
	fout[1].open("./Y1.txt", std::fstream::trunc);
	fout[2].open("./Y2.txt", std::fstream::trunc);
	fout[3].open("./Z.txt",  std::fstream::trunc);

	for(int i = 0; i < 4; ++i)
	{
		for(int k = 0; k <= Nt / sample_itv; ++k)
		{
			for(int j = 0; j <= Nx; ++j)
				fout[i] << uh_rec[i][k][j] << ' ';
			fout[i] << endl;
		}
		fout[i].close();
	}
	cout << "  Done." << endl;
	cout << "time elapsed: " << (double)(clock() - clock_begin) / CLOCKS_PER_SEC << "s" << endl;

	delete_1D(mesh);
	delete_1D(half_mesh);
	delete_1D(alpha1_half);
	delete_1D(alpha2_half);
	delete_1D(gamma1_half);
	delete_1D(gamma2_half);

	delete_2D(uh);
	delete_2D(uh_half);
	delete_3D(uh_rec);
	delete_2D(source);
	delete_2D(temp);

	return 0;
}
