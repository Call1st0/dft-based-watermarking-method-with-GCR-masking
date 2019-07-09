#pragma once

#include <armadillo>
#include <vector>
#include <string>
using namespace arma;
typedef void (*pfnPrintProgress)(size_t current, size_t total);

class Gamut
{
public:
	mat mat_Lab;
	Mat<int> mat_ch;
	size_t m_no_Lab;
	size_t m_no_ch;
	vec wpt;
	vec bpt;
	vec centroid;	

	Gamut(double *Lab, size_t no_Lab, int *ch, size_t no_ch);
	void map(double *Lab_set, double *Lab_set_new, size_t no_Lab_set, unsigned short *outgamarr, pfnPrintProgress PrintCallback);
	~Gamut();

private:
	std::vector<int> oct1;
	std::vector<int> oct2;
	std::vector<int> oct3;
	std::vector<int> oct4;
	std::vector<int> oct5;
	std::vector<int> oct6;
	std::vector<int> oct7;
	std::vector<int> oct8;
	std::map<int, std::vector<int>> vec_map;

	void PrintProgress(pfnPrintProgress PrintCallback, size_t current, size_t total);
	int GetOctand(vec Lab_pos);
};

