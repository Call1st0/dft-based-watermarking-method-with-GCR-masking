
#include"armadillo"
#include"Gamut.h"
#include"RevInt.h"
#include<string>

//Function declarations
extern "C"
{
	//Gamut mapping
	void *makegma(double *Lab, size_t no_Lab, int *ch, size_t no_ch);
	void map(void *gma, double *Lab_set, double *Lab_set_new, size_t no_Lab_set, unsigned short *outgamarr, pfnPrintProgress PrintCallback);
	void deletegma(void *gma);

	//Revese interpolation
	void *makeRevInt(const char* prof, int nogrdpts, int noclpts);
	void interp(void *rev, double *Lab, size_t no_Lab, double *kvals, size_t no_kvals);
	void interpex(void *rev, double *Lab, size_t no_Lab, double *kvals, size_t no_kvals);
	void interpLut(void *rev, double *lutcmykarr, unsigned short *outgamarr, double* kminmax, int nogrdpts);
	void getLabExtmRng(void *rev, double *labextmarr, double *rangearr);
	void *makegammapForm(float *LUTvals, unsigned int nogrdpts);
	void *makelab2kForm(float *LUTvals, unsigned int nogrdpts);
	void *makelabk2cmykForm(float *LUTvals, size_t nogrdpts);
	void applygammapForm(void *gammapform, double* LabIn, double* LabOut, size_t size);
	void applylab2kForm(void *lab2kform, double* LabIn, double* Kout, size_t size);
	void applylabk2cmykForm(void *labk2cmykform, double* LabKin, double* CMYKout, size_t size);
	void deletegammapForm(void *gammapform);
	void deletelab2kForm(void *labk2kform);
	void deletelabk2cmykForm(void *labk2cmykform);
	void deleteRevInt(void *rev);
}


//Function definitons
extern "C" void *makegma(double *Lab, size_t no_Lab, int *ch, size_t no_ch)
{
	Gamut *gma =  new Gamut(Lab, no_Lab, ch, no_ch);
	return reinterpret_cast<void *>(gma);
}

extern "C" void map(void *gma, double *Lab_set, double *Lab_set_new, size_t no_Lab_set, unsigned short *outgamarr, pfnPrintProgress PrintCallback)
{
	try
	{
		reinterpret_cast<Gamut*>(gma)->map(Lab_set, Lab_set_new, no_Lab_set, outgamarr, PrintCallback);
	}
	catch (const std::exception& ex)
	{
		std::ofstream out("GCR_py_wrap.txt");
		std::string str(ex.what());
		out << str;
		out.close();
	}	
}

extern "C" void deletegma(void *gma)
{
	reinterpret_cast<Gamut*>(gma)->~Gamut();
}

extern "C" void *makeRevInt(const char* prof, int nogrdpts, int noclpts)
{
	RevInt *rev = new RevInt(prof, nogrdpts, noclpts);
	return reinterpret_cast<void *>(rev);
}

extern "C" void interp(void *rev, double *Lab, size_t no_Lab, double *kvals, size_t no_kvals)
{
	vec kvalsvec(kvals, no_kvals);
	mat res(3, no_Lab*no_kvals);
	size_t count = 0;
	for (size_t i = 0; i < no_Lab * 3; i += 3)
	{
		vec labpt(Lab[i], 3);
		res(span(0, 2), span(count*no_kvals, count*no_kvals + no_kvals - 1)) = 
			reinterpret_cast<RevInt*>(rev)->interp(labpt, kvalsvec);
		count++;
	}
}

extern "C" void interpex(void *rev, double *Lab, size_t no_Lab, double *kvals, size_t no_kvals)
{
	vec kvalsvec(kvals, no_kvals);
	mat res(3, no_Lab*no_kvals);
	size_t count = 0;
	for (size_t i = 0; i < no_Lab * 3; i += 3)
	{
		vec labpt(Lab[i], 3);
		res(span(0, 2), span(count*no_kvals, count*no_kvals + no_kvals - 1)) =
			reinterpret_cast<RevInt*>(rev)->interpex(labpt, kvalsvec);
		count++;
	}
}

extern "C" void interpLut(void *rev, double *lutcmykarr, unsigned short *outgamarr, double* kminmax, int nogrdpts)
{
	reinterpret_cast<RevInt*>(rev)->interpLut(lutcmykarr, outgamarr, kminmax, nogrdpts);
}

extern "C" void getLabExtmRng(void *rev, double *labextmarr, double *rangearr)
{
	reinterpret_cast<RevInt*>(rev)->getLabExtmRng(labextmarr, rangearr);
}

extern "C" void *makegammapForm(float *LUTvals, unsigned int nogrdpts)
{
	return reinterpret_cast<void*>(RevInt::makegammapForm(LUTvals, nogrdpts));
}

extern "C" void *makelab2kForm(float *LUTvals, unsigned int nogrdpts)
{
	return reinterpret_cast<void*>(RevInt::makelab2kForm(LUTvals, nogrdpts));
}

extern "C" void *makelabk2cmykForm(float *LUTvals, size_t nogrdpts)
{
	return reinterpret_cast<void*>(RevInt::makelabk2cmykForm(LUTvals, nogrdpts));
}

extern "C" void applygammapForm(void *gammapform, double* LabIn, double* LabOut, size_t size)
{
	cmsDoTransform(reinterpret_cast<cmsHTRANSFORM>(gammapform), LabIn, LabOut, size);
}

extern "C" void applylab2kForm(void *lab2kform, double* LabIn, double* Kout, size_t size)
{
	cmsDoTransform(reinterpret_cast<cmsHTRANSFORM>(lab2kform), LabIn, Kout, size);
}

extern "C" void applylabk2cmykForm(void *labk2cmykform, double* LabKin, double* CMYKout, size_t size)
{
	cmsDoTransform(reinterpret_cast<cmsHTRANSFORM>(labk2cmykform), LabKin, CMYKout, size);
}

extern "C" void deletegammapForm(void *gammapform)
{
	cmsDeleteTransform(reinterpret_cast<cmsHTRANSFORM>(gammapform));
}

extern "C" void deletelab2kForm(void *labk2kform)
{
	cmsDeleteTransform(reinterpret_cast<cmsHTRANSFORM>(labk2kform));
}

extern "C" void deletelabk2cmykForm(void *labk2cmykform)
{
	cmsDeleteTransform(reinterpret_cast<cmsHTRANSFORM>(labk2cmykform));
}

extern "C" void deleteRevInt(void *rev)
{
	reinterpret_cast<RevInt*>(rev)->~RevInt();
}