#pragma once
#include <functional>
#include <vector>

class SplineInterpolator {	
	std::function<double(double)> f;
	int n;
	double a, b;
	std::vector<double> x_i, y_i, spline_y_i;
	std::vector<double> a_i, b_i, c_i, d_i;
	void FillXYVectors();
	void FillCVector(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d);
public:
	SplineInterpolator(const std::function<double(double)>& f, int n, const double& a, const double& b);
	void SetN(int n);
	void Interpolate();
	double FindError();
	void Plot();
};

