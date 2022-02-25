#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "Utils.h"

std::vector<double> x_i;
std::vector<double> y_i;
std::vector<double> lagrange_saved;

const int N_FOR_COMPARISON = 100000;

double FuncExact(double x) {
	return sin(x) + x*cos(x);
}

void FillVectorsEvenly(int N) {
	x_i.clear();
	y_i.clear();
	double h = (b - a) / N;
	for (int i = 0; i < N; ++i) {
		x_i.push_back(a + i * h);
		y_i.push_back(FuncExact(x_i[i]));
	}
}

void FillVectorsNotEvenly(int n0) {
	x_i.clear();
	y_i.clear();
	for (int i = 0; i < n0; ++i) {
		x_i.push_back((b + a) / 2 + cos((pi * ((double) 2 * i + 1)) / ((double) 2 * n0)) * (b - a) / 2);
		y_i.push_back(FuncExact(x_i[i]));
	}
}

double Lagrange(int n, double x) {
	double sum = 0;
	for (int i = 0; i <= n; ++i) {
		double mult = 1;
		for (int k = 0; k <= n; ++k) {
			if (i != k) {
				mult *= (x - x_i[k]) / (x_i[i] - x_i[k]);
			}
		}
		sum += y_i[i] * mult;
	}
	return sum;
}

double FindError(double n) {
	int N = N_FOR_COMPARISON;
	double max_error = 0;
	double h = (b - a) / N;
	for (int i = 0; i <= N; ++i) {
		double x_i = a + i * h;		
		double difference = abs(FuncExact(x_i) - Lagrange(n, x_i));
		if (difference > max_error) {
			max_error = difference;
		}
	}
	return max_error;
}

void Task1() {
	//1.1 - 1.3
	int min_error_n = 1;
	double min_error = N_FOR_COMPARISON;
	std::cout << "===== Lagrange Errors =====" << std::endl;
	std::ofstream error_vs_n("error_vs_n.txt");
	for (int n = 1; n <= 15; ++n) {
		FillVectorsEvenly(n + 1);
		double error = FindError(n);
		std::cout << "n = " << n << "\terror: " << std::fixed << std::setprecision(10) << error << std::endl;
		error_vs_n << n << " " << std::fixed << std::setprecision(10) << error << std::endl;
		if (error < min_error) {
			min_error_n = n;
			min_error = error;
		}
	}
	error_vs_n.close();
	std::cout << "===== Lagrange Errors =====" << std::endl;
	DrawGraphic("\"" + PATH + "error_vs_n.txt\" title \"Error vs n\" with lines");
	//1.4
	int N = N_FOR_COMPARISON;
	FillVectorsEvenly(min_error_n + 1);
	double error = FindError(min_error_n);
	double h = (b - a) / N;
	std::vector<double> lagrange;
	for (int i = 0; i <= N; ++i) {
		double x_i = a + i * h;
		lagrange.push_back(Lagrange(min_error_n, x_i));
	}
	FillVectorsEvenly(N + 1);
	std::ofstream exact_vs_lagrange_error("exact_vs_lagrange_error.txt");
	for (int i = 0; i < x_i.size(); ++i) {
		exact_vs_lagrange_error << x_i[i] << " " << abs(y_i[i] - lagrange[i]) << std::endl;
	}
	exact_vs_lagrange_error.close();
	//для 2
	std::ofstream exact_vs_lagrange("exact_vs_lagrange.txt");
	for (int i = 0; i < x_i.size(); ++i) {
		exact_vs_lagrange << x_i[i] << " " << y_i[i] << " " << lagrange[i] << std::endl;
	}
	exact_vs_lagrange.close();
	lagrange_saved = lagrange;
	DrawGraphic("\"" + PATH + "exact_vs_lagrange_error.txt\" using 1:2 title \"Exact vs Lagrange error\" with lines");
}

double Factorial(int n) {
	if (n == 0 || n == 1) {
		return 1;
	}
	return n * Factorial(n - 1);
}

void Task2() {
	//2.1 - 2.3
	std::cout << "===== Lagrange Not Evenly Grid Error (n0) =====" << std::endl;
	int n0 = 15;
	FillVectorsNotEvenly(n0 + 1);
	std::cout << "n0 = 15:" << std::endl;
	double real_error = FindError(n0);
	double theoretical_error = 518918400 * pow((b - a), n0 + 1) * pow(2, 1 - 2 * (n0+1)) / Factorial(n0 + 1);
	std::cout << "  Real error: " << std::fixed << std::setprecision(10) << real_error << std::endl;	
	std::cout << "  Theoretical error: " << std::fixed << std::setprecision(10) << theoretical_error << std::endl;
	std::cout << "  Difference (abs): " << std::fixed << std::setprecision(10) << abs(theoretical_error - real_error) << std::endl;	
	//2.4
	int N = N_FOR_COMPARISON;
	FillVectorsNotEvenly(n0 + 1);
	double h = (b - a) / N;
	double max_error = 0;
	std::ofstream lagrange_not_evenly_grid("lagrange_not_evenly_grid.txt");
	for (int i = 0; i <= N; ++i) {
		double x_i = a + i * h;
		lagrange_not_evenly_grid << x_i << " " << Lagrange(n0, x_i) << std::endl;
		double error = abs(lagrange_saved[i] - Lagrange(n0, x_i));
		if (error > max_error) {
			max_error = error;
		}
	}
	lagrange_not_evenly_grid.close();
	std::cout << "  Max error (grid comparison): " << std::fixed << std::setprecision(10) << max_error << std::endl;
	std::cout << "===== Lagrange Not Evenly Grid Error (n0) =====" << std::endl;
	DrawGraphic("\"" + PATH + "exact_vs_lagrange.txt\" using 1:3 title \"Evenly grid\" with lines, \"" + PATH + "lagrange_not_evenly_grid.txt\" using 1:2 title \"Not evenly grid\" with lines");
}

double DividedDifference(int start, int end) {
	if (start == end) {
		return y_i[start];
	}
	return (DividedDifference(start, end - 1) - DividedDifference(start + 1, end)) / (x_i[start] - x_i[end]);
}

double Newton(int n, double x) {
	double sum = y_i[0];
	for (int k = 1; k <= n; ++k) {
		double mult = 1;
		for (int i = 0; i <= k - 1; ++i) {
			mult *= x - x_i[i];
		}
		mult *= DividedDifference(0, k);
		sum += mult;
	}
	return sum;
}

void Task3() {
	//3.1 - 3.2
	std::cout << "===== Newton =====" << std::endl;
	int n0 = 15;
	int N = 100;
	FillVectorsNotEvenly(n0 + 1);
	double h = (b - a) / N;
	std::ofstream newton_evenly_grid("newton_evenly_grid.txt");
	double max_error = 0;
	for (int i = 0; i <= N; ++i) {
		double x_i = a + i * h;
		newton_evenly_grid << x_i << " " << Newton(n0, x_i) << std::endl;
		double error = abs(Newton(n0, x_i) - lagrange_saved[i]);
		if (error > max_error) {
			max_error = error;
		}
	}
	newton_evenly_grid.close();
	std::cout << "  Max error: " << max_error << std::endl;
	std::cout << "===== Newton =====" << std::endl;
	DrawGraphic("\"" + PATH + "exact_vs_lagrange.txt\" using 1:3 title \"Lagrange\" with lines, \"" + PATH + "newton_evenly_grid.txt\" using 1:2 title \"Newton\" with lines");
}
