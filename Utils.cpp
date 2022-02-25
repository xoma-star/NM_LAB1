#pragma once
#include <fstream>

void DrawGraphic(const std::string& data) {
	std::ofstream file("file");
	file << "plot " << data << "; pause mouse keypress" << "\n";
	file.close();
	std::system("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persist file");
}