#include "SplineInterpolator.h"
#include "Utils.h"
#include <iostream>
#include <fstream>

double f(double x) {
	return sin(x) + x*cos(x);
}

void Task6() {	
	int n = 0;
	double error = 10000;
	int max_n = 2000;
	SplineInterpolator spline_interpolator(f, n + 1, a, b);	
	std::ofstream spline_n_vs_error("spline_n_vs_error.txt");
	std::cout << "===== Spline interpolation, try to find n =====" << std::endl;
	while (error >= delta ) {
		++n;
		spline_interpolator.Interpolate();
		error = spline_interpolator.FindError();
		std::cout << "n = " << n << " error: " << error << std::endl;
		spline_n_vs_error << n << " " << error << std::endl;
		spline_interpolator.SetN(n);
		if (n >= max_n) {
			break;
		}
	}
	if (n < max_n) {
		std::cout << n << std::endl;
	}
	else {
		std::cout << "n not found" << std::endl;
	}
	std::cout << "===== Spline interpolation, try to find n =====" << std::endl;
	spline_n_vs_error.close();
	DrawGraphic("\"" + PATH + "spline_n_vs_error.txt\" using 1:2 title \"Spline, n vs error\" with lines");
	spline_interpolator.Plot();
}