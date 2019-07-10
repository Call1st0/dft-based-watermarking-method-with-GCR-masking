#pragma once
#include <armadillo>
#include <lcms2.h>

using namespace arma;

class RevInt
{
public:
	int m_nogrdpts;
	int m_noclpts;
	mat m_cmykgrid;
	mat m_labgrid;
	mat m_labextm;
	vec m_labrng;
	cmsHTRANSFORM m_cfAtoB;
	cmsHTRANSFORM m_cfBtoA;

	RevInt(const char* prof, int nogrdpts, int noclpts);
	RevInt(const char* prof, const char* cubelist);
	mat interp(vec labpt, vec kvals);
	mat interpex(vec labpt, vec kvals);
	void interpLut(double *lutcmykarr, unsigned short *outgamarr, double* kminmax, int nogrdpts);
	void getLabExtmRng(double *labextmarr, double *rangearr);
	static cmsHTRANSFORM makegammapForm(float *LUTvals, unsigned int nogrdpts);
	static cmsHTRANSFORM makelab2kForm(float *LUTvals, unsigned int nogrdpts);
	static cmsHTRANSFORM makelabk2cmykForm(float *LUTvals, size_t nogrdpts);
	static void applygammapForm(cmsHTRANSFORM gammapform, double* LabIn, double* LabOut, size_t size);
	static void applylab2kForm(cmsHTRANSFORM labk2kform, double* LabIn, double* Kout, size_t size);
	static void applylabk2cmykForm(cmsHTRANSFORM labk2cmykform, double* LabKin, double* CMYKout, size_t size);
	static void deletegammapForm(cmsHTRANSFORM gammapform);
	static void deletelab2kForm(cmsHTRANSFORM labk2kform);
	static void deletelabk2cmykForm(cmsHTRANSFORM labk2cmykform);
	~RevInt();

private:
	std::vector<Mat<int> > simplices;
	std::vector<std::vector<std::vector<int> > > m_cubelist;
	
	mat interpexLut(vec labpt, vec kvals, std::vector<arma::s32_mat> *sol_simp_vec,
		std::vector<std::vector<double> > *sol_k_vec, double *kminout, double *kmaxout);
	cmsHTRANSFORM makeCForm(const char *profname, const char *luttype);
	mat permRep(vec x, int nocols, const char* varies);
	Mat<int> permWoRep(int x[], int size);
	int factorial(int n);
	int lutInd4D(int i, int j, int k, int l, int nopts);
	int lutInd3D(int i, int j, int k, int nopts);
	Mat<int> intersect(Mat<int> A, Mat<int> B);
	vec overlap(vec c0, vec c1);
};

