#define _USE_MATH_DEFINES

/* cd "/Users/tianxuwang/Documents/Research/bold-shy pde/code" */
/* g++ -std=c++11 -o bold_shy_pde_2D bold_shy_pde_2D.cpp */
/* ./bold_shy_pde_2D */

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>

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

constexpr double a1 = -0.0034;
constexpr double b1 = 0.085;
constexpr double c1 = 4.5384;
constexpr double a2 = -0.0004;
constexpr double b2 = 0.0012;
constexpr double c2 = 0.0516;
constexpr double xi1 = 0.6;
constexpr double xi2 = 0.6;
constexpr double eta = 1;
constexpr double tau1 = 1.235;
constexpr double tau2 = 33.33;

constexpr double d11 = 0.00;
constexpr double d20 = 0.25;
constexpr double d30 = 0.003;
constexpr double d2 = 0.0159;
constexpr double d31 = 0.00004;
constexpr double beta1 = 0.000;
constexpr double beta2 = 0.00;
constexpr double T0 = 10;
constexpr double kappa = 0.1;

constexpr double A = 100;
constexpr double B = 15;
constexpr int    Nx = 200;
constexpr int    Ny = 200;
constexpr int    Nt = 50000;
// This is the actual number of grid points in time which is used in the solver.
// should be near or at least Nx^2
constexpr int    sample_itv = 100; 
// We record the solution for each *sample_itv* steps. 
// Nt in Matlab is the same as *Nt / sample_itv* in this program
constexpr double T_final = 500;
constexpr double dt = T_final / (double)Nt;
constexpr double dx = 2.0 * A / (double)Nx;
constexpr double dy = 2.0 * A / (double)Ny;
//
//static double T(double x, double y) { return T0 * exp(-(x * x + y * y) / (B * B)); }
//static double T(double x, double y) { return T0; }
//static double T(double x, double y) { return fmax(-(T0 / 100) * (x + 100) + T0,0); }
//static double T(double x, double y) { return fmax(T0 * (100 - x + y) / 300,0); }
static double T(double x, double y) { return fmax(fmax(T0 * exp(-(x * x + (y - 50) * (y - 50)) / (B * B)), T0 * exp(-((x + 50) * (x + 50) + (y + 50) * (y + 50)) / (B * B))), T0 * exp(-((x - 50) * (x - 50) + (y + 50) * (y + 50)) / (B * B))); }
//static double T(double x, double y) { return (T0 / (100*100*100*100)) * fmax(100*100*100*100 - (100 + x) * (100 - x) * (100 + y) * (100 - y),0); }
static double mu11(double T) { return fmax(a1 * T * T + b1 * T + c1, c1/2); }
static double mu21(double T) { return fmax(a2 * T * T + b2 * T + c2, c2/2); }
static double mu12(double T) { return xi1 * mu11(T); }
static double mu22(double T) { return xi2 * mu21(T); }
static double alpha1(double T) { return eta * mu11(T); }
static double alpha2(double T) { return eta * mu12(T); }
static double gamma1(double T) { return eta * mu21(T); }
static double gamma2(double T) { return eta * mu22(T); }
//static double d21(double T) { return d20; }
//static double d22(double T) { return d20; }
//static double d3(double T)  { return d30 + d31 * T; }
static double d21(double T) { return 0.013983 * (exp(0.160884 * T) - 1.0) + 0.25; }

static double d22(double T) { return xi1 * 0.013983 * (exp(0.160884 * T) - 1.0) + 0.25; }

static double d3(double T)  { return 0.010872 * (exp(0.071760 * T) - 1.0) + 0.003; }
static double h(double z)   { return z / (1 + kappa * z); }

static double X0 (double x, double y) { return 0.5 + 0.005 * cos(x * y / 1000 * 7); }
static double Y10(double x, double y) { return 0.5 - 0.005 * cos(x * y / 1000 * 5); }
static double Y20(double x, double y) { return 0.5 - 0.005 * cos(x * y / 1000 * 5); }
static double Z0 (double x, double y) { return 5 + 0.005 * cos(x * y / 1000 * 9); }

void interpolate_x_half(double ***uh, double ***uh_half)
{
	for(int ic = 0; ic < 4; ++ic)
	{
		for(int ix = 0; ix < Nx + 2; ++ix)
		for(int iy = 0; iy < Ny + 3; ++iy)
		{
			if(ix == 0 || ix == Nx + 1)
				uh_half[ic][ix][iy] = (uh[ic][ix][iy] + uh[ic][ix + 1][iy]) / 2.0;
			else
				uh_half[ic][ix][iy] = (9.0 * uh[ic][ix][iy] + 9.0 * uh[ic][ix + 1][iy] - uh[ic][ix - 1][iy] - uh[ic][ix + 2][iy]) / 16.0;
		}
	}
}

void interpolate_y_half(double ***uh, double ***uh_half)
{
	for(int ic = 0; ic < 4; ++ic)
	{
		for(int ix = 0; ix < Nx + 3; ++ix)
		for(int iy = 0; iy < Ny + 2; ++iy)
		{
			if(iy == 0 || iy == Ny + 1)
				uh_half[ic][ix][iy] = (uh[ic][ix][iy] + uh[ic][ix][iy + 1]) / 2.0;
			else
				uh_half[ic][ix][iy] = (9.0 * uh[ic][ix][iy] + 9.0 * uh[ic][ix][iy + 1] - uh[ic][ix][iy - 1] - uh[ic][ix][iy + 2]) / 16.0;
		}
	}
}

void setSourceVector(double ***uh, double *mesh_x, double *mesh_y, double ***source)
{
	for(int ix = 0; ix < Nx + 3; ++ix)
	for(int iy = 0; iy < Ny + 3; ++iy)
	{
		if(ix == 0 || ix == Nx + 2 || iy == 0 || iy == Ny + 2)
			for(int ic = 0; ic < 4; ++ic)
				source[ic][ix][iy] = 0;
		else
		{
			double Tval = T(mesh_x[ix], mesh_y[iy]);
			double fxy1 = mu11(Tval) * uh[0][ix][iy] * uh[1][ix][iy] / (1.0 + mu11(Tval) * tau1 * uh[0][ix][iy]);
			double fxy2 = mu12(Tval) * uh[0][ix][iy] * uh[2][ix][iy] / (1.0 + mu12(Tval) * tau1 * uh[0][ix][iy]);
			double gy1z = mu21(Tval) * uh[1][ix][iy] * uh[3][ix][iy] / (1.0 + mu21(Tval) * tau2 * uh[1][ix][iy]);
			double gy2z = mu22(Tval) * uh[2][ix][iy] * uh[3][ix][iy] / (1.0 + mu22(Tval) * tau2 * uh[2][ix][iy]);

			source[0][ix][iy] = r * uh[0][ix][iy] * (1 - uh[0][ix][iy] / K) - fxy1 - fxy2;
			source[1][ix][iy] = e1 * fxy1 - gy1z - d21(Tval) * uh[1][ix][iy] - beta1 * Tval * uh[1][ix][iy] + beta2 * Tval * uh[2][ix][iy];
			source[2][ix][iy] = e1 * fxy2 - gy2z - d22(Tval) * uh[2][ix][iy] + beta1 * Tval * uh[1][ix][iy] - beta2 * Tval * uh[2][ix][iy];
			source[3][ix][iy] = e2 * gy1z + e2 * gy2z - d3(Tval) * uh[3][ix][iy];
		}
	}
}

int main()
{
	clock_t clock_begin = clock();
	cout << "Solving PDE..." << std::endl;

	double  *mesh_x      = new_1D<double> (Nx + 3);
	double  *mesh_y      = new_1D<double> (Ny + 3);
	double  *half_mesh_x = new_1D<double> (Nx + 2);
	double  *half_mesh_y = new_1D<double> (Ny + 2);
	double **alpha1_x_half = new_2D<double> (Nx + 2, Ny + 3);
	double **alpha2_x_half = new_2D<double> (Nx + 2, Ny + 3);
	double **gamma1_x_half = new_2D<double> (Nx + 2, Ny + 3);
	double **gamma2_x_half = new_2D<double> (Nx + 2, Ny + 3);
	double **alpha1_y_half = new_2D<double> (Nx + 3, Ny + 2);
	double **alpha2_y_half = new_2D<double> (Nx + 3, Ny + 2);
	double **gamma1_y_half = new_2D<double> (Nx + 3, Ny + 2);
	double **gamma2_y_half = new_2D<double> (Nx + 3, Ny + 2);

	for(int ix = -1; ix <= Nx + 1; ++ix)
		mesh_x[ix + 1] = ix * dx - A;
	for(int iy = -1; iy <= Ny + 1; ++iy)
		mesh_y[iy + 1] = iy * dy - A;
	for(int ix = 0; ix <= Nx + 1; ++ix)
		half_mesh_x[ix] = (ix - 0.5) * dx - A;
	for(int iy = 0; iy <= Ny + 1; ++iy)
		half_mesh_y[iy] = (iy - 0.5) * dy - A;

	for(int ix = 0; ix <= Nx + 1; ++ix)
	for(int iy = 0; iy <= Ny + 2; ++iy)
	{
		double T_val = T(half_mesh_x[ix], mesh_y[iy]);
		alpha1_x_half[ix][iy] = alpha1(T_val);
		alpha2_x_half[ix][iy] = alpha2(T_val);
		gamma1_x_half[ix][iy] = gamma1(T_val);
        gamma2_x_half[ix][iy] = gamma2(T_val);
	}

	for(int ix = 0; ix <= Nx + 2; ++ix)
	for(int iy = 0; iy <= Ny + 1; ++iy)
	{
		double T_val = T(mesh_x[ix], half_mesh_y[iy]);
		alpha1_y_half[ix][iy] = alpha1(T_val);
		alpha2_y_half[ix][iy] = alpha2(T_val);
		gamma1_y_half[ix][iy] = gamma1(T_val);
        gamma2_y_half[ix][iy] = gamma2(T_val);
	}

	double ***uh        = new_3D<double> (4, Nx + 3, Ny + 3);
	double ***uh_x_half = new_3D<double> (4, Nx + 2, Ny + 3);
	double ***uh_y_half = new_3D<double> (4, Nx + 3, Ny + 2);
	double ***uh_rec = new_3D<double> (4, Nt / sample_itv + 1, (Nx + 1) * (Ny + 1));

	std::mt19937::result_type seed = 31415926; // time(0)
    auto real_rand = std::bind(std::uniform_real_distribution<double>(0., 1.), std::mt19937(seed));
	for(int i = 0; i < Nx + 3; ++i)
	for(int j = 0; j < Ny + 3; ++j)
	{
//		uh[0][i][j] = X0 (mesh_x[i], mesh_y[j]) + real_rand() * 0.001;
//		uh[1][i][j] = Y10(mesh_x[i], mesh_y[j]) + real_rand() * 0.001;
//		uh[2][i][j] = Y20(mesh_x[i], mesh_y[j]) + real_rand() * 0.001;
//		uh[3][i][j] = Z0 (mesh_x[i], mesh_y[j]) + real_rand() * 0.001;
		uh[0][i][j] = real_rand();
		uh[1][i][j] = real_rand();
		uh[2][i][j] = real_rand();
		uh[3][i][j] = real_rand();
	}

	for(int ic = 0; ic < 4; ++ic)
	{
		for(int ix = 0; ix <= Nx; ++ix)
		for(int iy = 0; iy <= Ny; ++iy)
			uh_rec[ic][0][ix * (Ny + 1) + iy] = uh[ic][ix + 1][iy + 1];
	}

	int iter = 0;
	double ***source = new_3D<double> (4, Nx + 3, Ny + 3);
	double ***temp   = new_3D<double> (4, Nx + 3, Ny + 3);

	while(iter < Nt)
	{
		interpolate_x_half(uh, uh_x_half);
		interpolate_y_half(uh, uh_y_half);
		setSourceVector(uh, mesh_x, mesh_y, source);

		for(int ix = 1; ix <= Nx + 1; ++ix)
		for(int iy = 1; iy <= Ny + 1; ++iy)
		{
			temp[0][ix][iy] = dt * source[0][ix][iy] + uh[0][ix][iy] 
				+ dt / dx / dx * delta1 * (uh[0][ix - 1][iy] + uh[0][ix + 1][iy] - 2.0 * uh[0][ix][iy])
				+ dt / dy / dy * delta1 * (uh[0][ix][iy - 1] + uh[0][ix][iy + 1] - 2.0 * uh[0][ix][iy]);
			temp[1][ix][iy] = dt * source[1][ix][iy] + uh[1][ix][iy] 
				+ dt / dx / dx * (delta2 * (uh[1][ix - 1][iy] + uh[1][ix + 1][iy] - 2.0 * uh[1][ix][iy])
					+ alpha1_x_half[ix - 1][iy] * h(uh_x_half[1][ix - 1][iy]) * (uh[0][ix][iy] - uh[0][ix - 1][iy])
					+ alpha1_x_half[ix][iy]     * h(uh_x_half[1][ix][iy])     * (uh[0][ix][iy] - uh[0][ix + 1][iy]) )
				+ dt / dy / dy * (delta2 * (uh[1][ix][iy - 1] + uh[1][ix][iy + 1] - 2.0 * uh[1][ix][iy])
					+ alpha1_y_half[ix][iy - 1] * h(uh_y_half[1][ix][iy - 1]) * (uh[0][ix][iy] - uh[0][ix][iy - 1])
					+ alpha1_y_half[ix][iy]     * h(uh_y_half[1][ix][iy])     * (uh[0][ix][iy] - uh[0][ix][iy + 1]) );
			temp[2][ix][iy] = dt * source[2][ix][iy] + uh[2][ix][iy]
				+ dt / dx / dx * (delta2 * (uh[2][ix - 1][iy] + uh[2][ix + 1][iy] - 2.0 * uh[2][ix][iy])
					+ alpha2_x_half[ix - 1][iy] * h(uh_x_half[2][ix - 1][iy]) * (uh[0][ix][iy] - uh[0][ix - 1][iy])
					+ alpha2_x_half[ix][iy]     * h(uh_x_half[2][ix][iy])     * (uh[0][ix][iy] - uh[0][ix + 1][iy]) )
				+ dt / dy / dy * (delta2 * (uh[2][ix][iy - 1] + uh[2][ix][iy + 1] - 2.0 * uh[2][ix][iy])
					+ alpha2_y_half[ix][iy - 1] * h(uh_y_half[2][ix][iy - 1]) * (uh[0][ix][iy] - uh[0][ix][iy - 1])
					+ alpha2_y_half[ix][iy]     * h(uh_y_half[2][ix][iy])     * (uh[0][ix][iy] - uh[0][ix][iy + 1]));
			temp[3][ix][iy] = dt * source[3][ix][iy] + uh[3][ix][iy]
				+ dt / dx / dx * ( delta3 * (uh[3][ix - 1][iy] + uh[3][ix + 1][iy] - 2.0 * uh[3][ix][iy])
					+ gamma1_x_half[ix - 1][iy] * h(uh_x_half[3][ix - 1][iy]) * (uh[1][ix][iy] - uh[1][ix - 1][iy])
					+ gamma1_x_half[ix][iy]     * h(uh_x_half[3][ix][iy])     * (uh[1][ix][iy] - uh[1][ix + 1][iy])
					+ gamma2_x_half[ix - 1][iy] * h(uh_x_half[3][ix - 1][iy]) * (uh[2][ix][iy] - uh[2][ix - 1][iy])
					+ gamma2_x_half[ix][iy]     * h(uh_x_half[3][ix][iy])     * (uh[2][ix][iy] - uh[2][ix + 1][iy]) )
				+ dt / dy / dy * ( delta3 * (uh[3][ix][iy - 1] + uh[3][ix][iy + 1] - 2.0 * uh[3][ix][iy])
					+ gamma1_y_half[ix][iy - 1] * h(uh_y_half[3][ix][iy - 1]) * (uh[1][ix][iy] - uh[1][ix][iy - 1])
					+ gamma1_y_half[ix][iy]     * h(uh_y_half[3][ix][iy])     * (uh[1][ix][iy] - uh[1][ix][iy + 1])
					+ gamma2_y_half[ix][iy - 1] * h(uh_y_half[3][ix][iy - 1]) * (uh[2][ix][iy] - uh[2][ix][iy - 1])
					+ gamma2_y_half[ix][iy]     * h(uh_y_half[3][ix][iy])     * (uh[2][ix][iy] - uh[2][ix][iy + 1]) );
            temp[0][ix][iy] = fmax(temp[0][ix][iy], 0.);
            temp[1][ix][iy] = fmax(temp[1][ix][iy], 0.);
            temp[2][ix][iy] = fmax(temp[2][ix][iy], 0.);
            temp[3][ix][iy] = fmax(temp[3][ix][iy], 0.);
		}

		memcpy(uh[0][0], temp[0][0], (size_t)sizeof(double) * 4 * (Nx + 3) * (Ny + 3));
		for(int ic = 0; ic < 4; ++ic)
		{
			for(int iy = 0; iy <= Ny + 2; ++iy)
			{
				uh[ic][0][iy] = uh[ic][2][iy];
				uh[ic][Nx + 2][iy] = uh[ic][Nx][iy];
			}
			for(int ix = 0; ix <= Nx + 2; ++ix)
			{
				uh[ic][ix][0] = uh[ic][ix][2];
				uh[ic][ix][Ny + 2] = uh[ic][ix][Ny];
			}
		}

		++iter;
		if(iter % sample_itv == 0)
		{
			for(int ic = 0; ic < 4; ++ic)
				for(int ix = 0; ix <= Nx; ++ix)
				for(int iy = 0; iy <= Ny; ++iy)
					uh_rec[ic][iter / sample_itv][ix * (Ny + 1) + iy] = uh[ic][ix + 1][iy + 1];
//			cout << "  t = " << dt * iter << " / " << T_final << endl;
		}
	}
//	cout << "Done." << endl;
//	cout << "Time elapsed: " << (double)(clock() - clock_begin) / CLOCKS_PER_SEC << "s" << endl;

	cout << "Saving data to file..." << std::flush;
	std::ofstream fout[4];
	fout[0].open("./X.txt",  std::fstream::trunc);
	fout[1].open("./Y1.txt", std::fstream::trunc);
	fout[2].open("./Y2.txt", std::fstream::trunc);
	fout[3].open("./Z.txt",  std::fstream::trunc);

	for(int ic = 0; ic < 4; ++ic)
	{
		for(int it = 0; it <= Nt / sample_itv; ++it)
		{
			for(int ix = 0; ix <= Nx; ++ix)
			for(int iy = 0; iy <= Ny; ++iy)
				fout[ic] << uh_rec[ic][it][ix * (Ny + 1) + iy] << ' ';
			fout[ic] << endl;
		}
		fout[ic].close();
	}
	cout << "  Done." << endl;
	cout << "time elapsed: " << (double)(clock() - clock_begin) / CLOCKS_PER_SEC << "s" << endl;

	delete_1D(mesh_x);
	delete_1D(mesh_y);
	delete_1D(half_mesh_x);
	delete_1D(half_mesh_y);

	delete_2D(alpha1_x_half);
	delete_2D(alpha2_x_half);
	delete_2D(gamma1_x_half);
	delete_2D(gamma2_x_half);
	delete_2D(alpha1_y_half);
	delete_2D(alpha2_y_half);
	delete_2D(gamma1_y_half);
	delete_2D(gamma2_y_half);

	delete_3D(uh);
	delete_3D(uh_x_half);
	delete_3D(uh_y_half);
	delete_3D(uh_rec);
	delete_3D(source);
	delete_3D(temp);

	return 0;
}
