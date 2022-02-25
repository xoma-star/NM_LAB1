#include <vector>
#include <cmath>
#include <vector>
#include <fstream>

std::vector<double> t_i;
std::vector<double> a_k;
std::vector<double> b_k;

double g(double t) {
	return sin(t) + t*cos(t);
}

double Replacement(double point) {
	//[0, 3] -> [0, 2pi]
	return (((point - a) / (b - a)) * 2 * pi);
	//return point;
}

double ReverseReplacement(double point) {
	//[0, 2pi] -> [0, 3]
	return ((b * point / (2 * pi)) - a * (point / (2 * pi) - 1));
	//return point;
}

void FillVectors(int n) {
	t_i.clear();
	for (int i = 1; i <= 2 * n + 1; ++i) {		
		t_i.push_back((2 * pi * (i - 1)) / (2 * n + 1));
	}
}

void CalculateCoefficients(int n) {
	a_k.clear();
	b_k.clear();
	a_k.resize(n + 1, 0);
	b_k.resize(n + 1, 0);
	for (int i = 1; i <= 2 * n + 1; ++i) {
		double x_i = t_i[i - 1];
		double y_i = g(ReverseReplacement(x_i));
		a_k[0] += y_i;		
		for (int k = 1; k <= n; ++k) {
			a_k[k] += y_i * cos(k * x_i);
			b_k[k] += y_i * sin(k * x_i);
			if (i == 2 * n + 1) {
				a_k[k] *= 2. / (2 * n + 1);
				b_k[k] *= 2. / (2 * n + 1);
			}
		}
	}
	a_k[0] *= 2. / (2 * n + 1);
}

double TrigonometricInterpolation(int n, double t) {
	CalculateCoefficients(n);
	t = Replacement(t);
	double sum = a_k[0] / 2;
	for (int k = 1; k <= n; ++k) {
		sum += a_k[k] * cos(k * t) + b_k[k] * sin(k * t);
	}
	return sum;
}

double FindErrorTrigonometric(double n) {
	int N = 100;
	double max_error = 0;
	double h = (b - a) / N;
	for (int i = 0; i <= N; ++i) {
		double x_i = a + i * h;
		double difference = abs(g(x_i) - TrigonometricInterpolation(n, x_i));
		if (difference > max_error) {
			max_error = difference;
		}
	}
	return max_error;
}

void Plot() {
	int N = 100;
	int n0 = 105;
	FillVectors(n0 + 1);
	double h = (b - a) / N;
	std::ofstream trigonometric("trigonometric.txt");
	for (int i = 0; i <= N; ++i) {
		double x_i = a + i * h;
		trigonometric << x_i << " " << g(x_i) << " " << TrigonometricInterpolation(n0, x_i) << std::endl;
	}
	trigonometric.close();
	DrawGraphic("\"" + PATH + "trigonometric.txt\" using 1:2 title \"Exact\" with lines, \"" + PATH + "trigonometric.txt\" using 1:3 title \"Trigonometric\" with lines");
}

void Task4() {
	int n = 1;
	double error = 10000;
	int max_n = 500;
	std::ofstream trigonometric_n_vs_error("trigonometric_n_vs_error.txt");	
	std::cout << "===== Trigonometric interpolation, try to find n =====" << std::endl;
	while (error >= new_delta) {
		FillVectors(n);
		error = FindErrorTrigonometric(n);
		std::cout << "n = " << n << " error: " << error << std::endl;		
		trigonometric_n_vs_error << n << " " << error << std::endl;		
		if (++n >= max_n) {
			break;
		}
	}
	if (n < max_n) {
		std::cout << n << std::endl;
	}
	else {
		std::cout << "n not found" << std::endl;
	}
	std::cout << "===== Trigonometric interpolation, try to find n =====" << std::endl;
	trigonometric_n_vs_error.close();
	DrawGraphic("\"" + PATH + "trigonometric_n_vs_error.txt\" using 1:2 title \"Trigonometric, n vs error\" with lines");
	Plot();
}
