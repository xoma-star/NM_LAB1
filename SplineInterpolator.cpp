#include "SplineInterpolator.h"
#include <fstream>
#include "Utils.h"
#include <iostream>

void SplineInterpolator::FillXYVectors() {
	x_i.clear();
	y_i.clear();
	double h = (b - a) / n;
	for (int i = 0; i <= n; ++i) {
		x_i.push_back(a + i * h);
		y_i.push_back(f(x_i[i]));
	}
}

void SplineInterpolator::FillCVector(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) {
	//tridiagonal method
	std::vector<double> dzeta(n + 2, 0);
	std::vector<double> eta(n + 2, 0);
	for (int i = 1; i <= n; ++i) {
		if (b[i] - a[i] * dzeta[i] != 0) {
			dzeta[i + 1] = c[i] / (b[i] - a[i] * dzeta[i]);
			eta[i + 1] = (a[i] * eta[i] - d[i]) / (b[i] - a[i] * dzeta[i]);
		}
	}
	for (int i = n; i >= 1; --i) {
		double value = dzeta[i + 1] * c_i[i + 1] + eta[i + 1];
		c_i[i] = dzeta[i + 1] * c_i[i + 1] + eta[i + 1];
	}
}

SplineInterpolator::SplineInterpolator(const std::function<double(double)>& f, int n, const double& a, const double& b) : f(f), n(n), a(a), b(b) {

}

void SplineInterpolator::SetN(int n) {
	this->n = n;
}

void SplineInterpolator::Interpolate() {
	FillXYVectors();
	a_i.clear();
	b_i.clear();
	c_i.clear();
	d_i.clear();
	a_i.resize(n + 2, 0);
	b_i.resize(n + 2, 0);
	c_i.resize(n + 2, 0);
	d_i.resize(n + 2, 0);
	std::vector<double> a(n + 2, 0);
	std::vector<double> b(n + 2, 0);
	std::vector<double> c(n + 2, 0);
	std::vector<double> d(n + 2, 0);
	for (int i = 2; i <= n; ++i) {
		double h_i_prev = x_i[i - 1] - x_i[i - 2];
		double h_i = x_i[i] - x_i[i - 1];
		a[i] = h_i_prev;
		b[i] = 2 * (h_i_prev + h_i);
		c[i] = h_i;
		d[i] = 3 * ((y_i[i] - y_i[i - 1]) / h_i - (y_i[i - 1] - y_i[i - 2]) / h_i_prev);
	}
	FillCVector(a, b, c, d);
	for (int i = 1; i <= n - 1; ++i) {
		double h_i = x_i[i] - x_i[i - 1];
		b_i[i] = (y_i[i] - y_i[i - 1]) / h_i - 1. / 3 * h_i * (c_i[i + 1] + 2 * c_i[i]);
		d_i[i] = (c_i[i + 1] - c_i[i]) / (3 * h_i);
		a_i[i] = y_i[i - 1];
	}
	double h_n = x_i[n] - x_i[n - 1];
	b_i[n] = (y_i[n] - y_i[n - 1]) / h_n - 2. / 3 * h_n * c_i[n];
	d_i[n] = -c_i[n] / (3 * h_n);
	a_i[n] = y_i[n - 1];
}

double SplineInterpolator::FindError() {
	double max_error = abs(f(a) - a_i[1]);
	double h = (b - a) / n;
	for (int i = 1; i <= n; ++i) {
		double x_i_prev = a + (i - 1) * h;
		for (int j = 0; j <= 10; ++j) {
			double x_j = a + (i - 1) * h + (h / 10) * j;
			double phi = a_i[i] + b_i[i] * (x_j - x_i_prev) + c_i[i] * pow((x_j - x_i_prev), 2) + d_i[i] * pow((x_j - x_i_prev), 3);
			double difference = abs(f(x_j) - phi);
			if (difference > max_error) {
				max_error = difference;
			}
		}
	}
	return max_error;
}

void SplineInterpolator::Plot() {
	double h = (b - a) / n;
	std::ofstream spline("spline.txt");
	spline << a << " " << a_i[1] << std::endl;
	for (int i = 1; i <= n; ++i) {		
		double x_i_prev = a + (i - 1) * h;
		for (int j = 0; j <= 10; ++j) {
			double x_j = a + (i - 1) * h + (h / 10) * j;
			double phi = a_i[i] + b_i[i] * (x_j - x_i_prev) + c_i[i] * pow((x_j - x_i_prev), 2) + d_i[i] * pow((x_j - x_i_prev), 3);
			spline << x_j << " " << phi << std::endl;
		}		
	}
	spline.close();
	DrawGraphic("\"" + PATH + "lagrange_not_evenly_grid.txt\" using 1:2 title \"Exact\" with lines, \"" + PATH + "spline.txt\" using 1:2 title \"Spline\" with lines");
}
