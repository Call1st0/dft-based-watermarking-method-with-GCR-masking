#include "lcmswrap.h"
#include <string>
#include <unistd.h>

//Function declarations
extern "C"
{
	void *makecform(const char *profname, const char *luttype);
	void applycform(void *cform, double *input, double *output, size_t length);
	void deletecform(void *cform);
	void *makeabscform(float *LUTvals, unsigned int nogrdpts);
	void applyabscform(void *abscform, double* LabKin, double* CMYKout, size_t size);
	void deleteabscform(void *abscform);
	int testadd(int a, int b);
}

//Function definitons


extern "C" void *makecform(const char *profname, const char *luttype)
{
	cmsHTRANSFORM hTransform = NULL;
	cmsHPROFILE hInProfile;
	cmsHPROFILE hOutProfile;
	hInProfile = cmsOpenProfileFromFile(profname, "r");
	int counter = 0;
	while (!hInProfile && counter < 5)
	{
		usleep(1000);
		hInProfile = cmsOpenProfileFromFile(profname, "r");
		counter = counter + 1;
	}
	if (!hInProfile) {
		return NULL;
	}
	if (cmsGetDeviceClass(hInProfile) != 0x70727472) {
		return NULL;
	}
	hOutProfile = cmsCreateLab4Profile(NULL);

	std::string luttype_s(luttype);
	if (!luttype_s.compare("AtoB1")) {
		hTransform = cmsCreateTransform(hInProfile, TYPE_CMYK_DBL, hOutProfile, TYPE_Lab_DBL, INTENT_RELATIVE_COLORIMETRIC, 0);
	}
	else if (!luttype_s.compare("BtoA1")) {
		hTransform = cmsCreateTransform(hOutProfile, TYPE_Lab_DBL, hInProfile, TYPE_CMYK_DBL, INTENT_RELATIVE_COLORIMETRIC, 0);
	}

	cmsCloseProfile(hInProfile);
	cmsCloseProfile(hOutProfile);

	return reinterpret_cast<void *>(hTransform);	
}

extern "C" void applycform(void *cform, double *input, double *output, size_t length)
{
	//Stara linija, daje CIE Lab u cudnom rasponu...
	cmsDoTransform(reinterpret_cast<cmsHTRANSFORM>(cform), input, output, length);

	//Umjesto toga, sljedece rjesava jednu po jednu koristeci cmsCIELab Lab varijablu
	//Ili ipak ne, greska je bila u Pythonu. Pointer na niz mora biti osvjezen
	//neposredno prije prosljedjivanja.
	/*cmsCIELab Lab;
	size_t LabCounter = 0;
	for (size_t i = 0; i < length * 4; i = i + 4)
	{
		cmsDoTransform(reinterpret_cast<cmsHTRANSFORM>(cform), input + i, &Lab, 1);
		output[LabCounter] = Lab.L;
		output[LabCounter + 1] = Lab.a;
		output[LabCounter + 2] = Lab.b;
		LabCounter = LabCounter + 3;
	}*/
}

extern "C" void deletecform(void *cform)
{
	cmsDeleteTransform(reinterpret_cast<cmsHTRANSFORM>(cform));
}

extern "C" void *makeabscform(float *LUTvals, unsigned int nogrdpts)
{
	cmsHTRANSFORM hTransformp;
	cmsHPROFILE hProfilep = cmsCreateProfilePlaceholder(0);

	cmsSetProfileVersion(hProfilep, 4.2);
	cmsSetDeviceClass(hProfilep, cmsSigAbstractClass);
	cmsSetColorSpace(hProfilep, cmsSigCmykData);
	cmsSetPCS(hProfilep, cmsSigCmykData);

	cmsFloat32Number* cmsf = LUTvals;
	cmsStage* cmsStageAllocCLutFloat(cmsContext ContextID,
		cmsUInt32Number nGridPoints,
		cmsUInt32Number inputChan,
		cmsUInt32Number outputChan,
		const cmsFloat32Number * Table);
	cmsPipeline* Pipeline = cmsPipelineAlloc(0, 4, 4);
	cmsPipelineInsertStage(Pipeline, cmsAT_BEGIN, cmsStageAllocCLutFloat(0, nogrdpts, 4, 4, cmsf));
	cmsWriteTag(hProfilep, cmsSigAToB1Tag, Pipeline);
	cmsPipelineFree(Pipeline);

	hTransformp = cmsCreateMultiprofileTransform(&hProfilep,
		1,
		TYPE_CMYK_DBL,
		TYPE_CMYK_DBL,
		INTENT_RELATIVE_COLORIMETRIC,
		0);
	cmsCloseProfile(hProfilep);

	return reinterpret_cast<void *>(hTransformp);
}

extern "C" void applyabscform(void *abscform, double* LabKin, double* CMYKout, size_t size)
{
	for (size_t i = 0; i < size * 4; ++i)
	{
		LabKin[i] = LabKin[i] * 100;
	}
	cmsDoTransform(reinterpret_cast<cmsHTRANSFORM>(abscform), LabKin, CMYKout, size);
	for (size_t i = 0; i < size * 4; ++i)
	{
		CMYKout[i] = CMYKout[i] / 100;
	}
}

extern "C" void deleteabscform(void *abscform)
{
	cmsDeleteTransform(reinterpret_cast<cmsHTRANSFORM>(abscform));
}

extern "C" int testadd(int a, int b)
{
	return int(a+b);
}