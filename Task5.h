#pragma once
#include "Utils.h"
#include <cmath>
#include <vector>
#include <iostream>

std::vector<double> x;
std::vector<double> prev_x;
std::vector<double> y;
std::vector<double> prev_y;
std::vector<double> lagrange;
std::vector<double> prev_lagrange;

double Taylor(double x) {
	double value = 0.1053992246 - 0.3161976738 * (x - 1.5) + 0.3688972861 * pow((x - 1.5), 2) - 0.1580988369 * pow((x - 1.5), 3)
		- 0.06587451538 * pow((x - 1.5), 4) + 0.1027642440 * pow((x - 1.5), 5) - 0.02942395020 * pow((x - 1.5), 6)
		- 0.01675094819 * pow((x - 1.5), 7) + 0.01363759313 * pow((x - 1.5), 8) - 0.0008234314422 * pow((x - 1.5), 9)
		- 0.002480489192 * pow((x - 1.5), 10) + 0.0008262118600 * pow((x - 1.5), 11) + 0.0002068619003 * pow((x - 1.5), 12)
		- 0.0001748468786 * pow((x - 1.5), 13) + 7.915488213 * pow(10, -6) * pow((x - 1.5), 14) + 0.00002172981950 * pow((x - 1.5), 15)
		- 5.063777182 * pow(10, -6) * pow((x - 1.5), 16) - 1.662841614 * pow(10, -6) * pow((x - 1.5), 17) + 8.397821782 * pow(10, -7) * pow((x - 1.5), 18)
		+ 4.243877339 * pow(10, -8) * pow((x - 1.5), 19) - 9.034403383 * pow(10, -8) * pow((x - 1.5), 20);
	return value;
}

double LagrangeInterpolation(int n, double point) {
	double sum = 0;
	for (int i = 0; i <= n; ++i) {
		double mult = 1;
		for (int k = 0; k <= n; ++k) {
			if (i != k) {				
				mult *= (point - prev_x[k]) / (prev_x[i] - prev_x[k]);				
			}
		}
		sum += prev_y[i] * mult;
	}
	return sum;
}

void Task5() {	
	std::cout << "===== Best Uniform Approximation Polynomial =====" << std::endl;
	double err = 0;
	int N = 100;
	int n = 19;
	double h = (b - a) / N;
	double max_taylor_err = 0;
	for (int i = 0; i < N; ++i) {
		double x_i = a + i * h;
		double difference = abs(Taylor(x_i) - exp(-pow(x_i, 2)));
		if (difference >= max_taylor_err) {
			max_taylor_err = difference;
		}
	}
	std::ofstream buap_n_vs_err("buap_n_vs_err.txt");
	std::cout << "n = 20 error = " << std::fixed << std::setprecision(10) << max_taylor_err << std::endl;
	buap_n_vs_err << 20 << " " << max_taylor_err << std::endl;
	//Векторы для построения ИМЛ для разложения Тейлора (20 точек для 19 степени)
	for (int i = 0; i <= n; ++i) {
		prev_x.push_back((b + a) / 2 + cos((pi * ((double)2 * i + 1)) / ((double)2 * (n + 1))) * (b - a) / 2);
		prev_y.push_back(Taylor(prev_x[i]));
	}
	double max_19_err = 0;
	//Построение ИМЛ для разложения Тейлора (19 степень)
	for (int i = 0; i < N; ++i) {
		double x_i = a + i * h;
		double value = LagrangeInterpolation(n, x_i);
		double difference = abs(Taylor(x_i) - value);
		if (difference >= max_19_err) {
			max_19_err = difference;
		}
		lagrange.push_back(value);
		//Векторы для построения ИМЛ для ИМЛ для разложения Тейлора (т.е. 18 степень, 19 точек, строим по ИМЛ 19 степени)
		if (i <= n - 1) { //i = 0 .. 19 - 1 -> 19 точек			
			x.push_back((b + a) / 2 + cos((pi * ((double)2 * i + 1)) / ((double)2 * (n + 1))) * (b - a) / 2);			
			y.push_back(LagrangeInterpolation(n, x[i]));
		}
	}
	std::cout << "n = 19 error = " << std::fixed << std::setprecision(10) << max_19_err << std::endl;
	buap_n_vs_err << 19 << " " << max_19_err << std::endl;
	while (err < delta) {
		n--; //18 .. 
		//Меняем местами векторы
		prev_x = x;
		x.clear();	
		prev_y = y;
		y.clear();
		prev_lagrange = lagrange;
		lagrange.clear();
		double max_error = 0;
		for (int i = 0; i < N; ++i) {
			double x_i = a + i * h;		
			//Построение ИМЛ 18... степени
			lagrange.push_back(LagrangeInterpolation(n, x_i));
			if (i <= n - 1) { //0..17 -> 18 точек, строим по ИМЛ 18 степени
				x.push_back((b + a) / 2 + cos((pi * ((double)2 * i + 1)) / ((double)2 * (n + 1))) * (b - a) / 2);
				y.push_back(LagrangeInterpolation(n, x[i]));
			}
			double difference = abs(prev_lagrange[i] - lagrange[i]);
			if (difference >= max_error) {
				max_error = difference;
			}
		}
		err = max_error;
		std::cout << "n = " << n << " error = " << std::fixed << std::setprecision(10) << err << std::endl;
		buap_n_vs_err << n << " " << err << std::endl;
	}
	n--;
	std::cout << "Optimal n: " << n << std::endl;
	buap_n_vs_err.close();
	DrawGraphic("\"" + PATH + "buap_n_vs_err.txt\" title \"BUAP Error vs n\" with lines");
	std::cout << "===== Best Uniform Approximation Polynomial =====" << std::endl;
}
