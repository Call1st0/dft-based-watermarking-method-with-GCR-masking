#include "RevInt.h"


RevInt::RevInt(const char* prof, int nogrdpts, int noclpts)
{
	m_nogrdpts = nogrdpts;
	m_noclpts = noclpts;
	m_cfAtoB = makeCForm(prof, "AtoB1");
	m_cfBtoA = makeCForm(prof, "BtoA1");
	m_cmykgrid = permRep(arma::linspace(0, 100, nogrdpts), 4, "first");
	m_cmykgrid = m_cmykgrid.t();
	double *labgrid = new double[m_cmykgrid.n_elem/4*3];
	cmsDoTransform(m_cfAtoB, &m_cmykgrid(0,0), labgrid, m_cmykgrid.n_elem/4);
	m_labgrid = mat(labgrid, 3, m_cmykgrid.n_elem/4);
	delete[] labgrid;

	//4D cube vertices
	Mat<int> cv = arma::conv_to<arma::s32_mat>::from(permRep(arma::linspace(0, 1, 2), 4, "first"));

	//4D simplices around main diagonal, 4!=24
	int x[] = { 0, 1, 2, 3 };
	Mat<int> conds = permWoRep(x, 4);
	simplices.resize(24);
	for (int i = 0; i < conds.n_rows; i++)
	{
		Mat<int> set1 = cv.rows(arma::find(cv.col(conds(i, 0)) >= cv.col(conds(i, 1)) == 1));
		Mat<int> set2 = cv.rows(arma::find(cv.col(conds(i, 1)) >= cv.col(conds(i, 2)) == 1));
		Mat<int> set3 = cv.rows(arma::find(cv.col(conds(i, 2)) >= cv.col(conds(i, 3)) == 1));
		Mat<int> set4 = cv.rows(arma::find(cv.col(conds(i, 0)) >= cv.col(conds(i, 2)) == 1));
		Mat<int> set5 = cv.rows(arma::find(cv.col(conds(i, 0)) >= cv.col(conds(i, 3)) == 1));
		Mat<int> set6 = cv.rows(arma::find(cv.col(conds(i, 1)) >= cv.col(conds(i, 3)) == 1));
		Mat<int> set = intersect(intersect(intersect(intersect(intersect(set1, set2), set3), set4), set5), set6);
		simplices[i] = set;
	}
	
	//Get L*a*b* min and max from grid
	m_labextm.resize(3,2);
	m_labextm.col(0) = arma::min(m_labgrid, 1);
	m_labextm.col(1) = arma::max(m_labgrid, 1);
	m_labrng = m_labextm.col(1) - m_labextm.col(0);

	//Initialize list of L*a*b* cubes which correspond to 4D CMYK cubes
	m_cubelist.resize(pow(noclpts-1, 3));

	//Place each cube in 3D L*a*b* cube list
	mat cubelab(3, 16);
	mat labextmcur(3, 2); //L*a*b* min and max for current cube
	for (int l = 0; l < nogrdpts - 1; l++)
	{
		for (int k = 0; k < nogrdpts - 1; k++)
		{
			for (int j = 0; j < nogrdpts - 1; j++)
			{
				for (int i = 0; i < nogrdpts - 1; i++)
				{
					cubelab.col(0) = m_labgrid.col(lutInd4D(i, j, k, l, nogrdpts));
					cubelab.col(1) = m_labgrid.col(lutInd4D(i+1, j, k, l, nogrdpts));
					cubelab.col(2) = m_labgrid.col(lutInd4D(i, j+1, k, l, nogrdpts));
					cubelab.col(3) = m_labgrid.col(lutInd4D(i+1, j+1, k, l, nogrdpts));
					cubelab.col(4) = m_labgrid.col(lutInd4D(i, j, k+1, l, nogrdpts));
					cubelab.col(5) = m_labgrid.col(lutInd4D(i+1, j, k+1, l, nogrdpts));
					cubelab.col(6) = m_labgrid.col(lutInd4D(i, j+1, k+1, l, nogrdpts));
					cubelab.col(7) = m_labgrid.col(lutInd4D(i+1, j+1, k+1, l, nogrdpts));
					cubelab.col(8) = m_labgrid.col(lutInd4D(i, j, k, l+1, nogrdpts));
					cubelab.col(9) = m_labgrid.col(lutInd4D(i+1, j, k, l+1, nogrdpts));
					cubelab.col(10) = m_labgrid.col(lutInd4D(i, j+1, k, l+1, nogrdpts));
					cubelab.col(11) = m_labgrid.col(lutInd4D(i+1, j+1, k, l+1, nogrdpts));
					cubelab.col(12) = m_labgrid.col(lutInd4D(i, j, k+1, l+1, nogrdpts));
					cubelab.col(13) = m_labgrid.col(lutInd4D(i+1, j, k+1, l+1, nogrdpts));
					cubelab.col(14) = m_labgrid.col(lutInd4D(i, j+1, k+1, l+1, nogrdpts));
					cubelab.col(15) = m_labgrid.col(lutInd4D(i+1, j+1, k+1, l+1, nogrdpts));
					//Get L*a*b* min and max for current cube
					labextmcur.col(0) = arma::min(cubelab, 1);
					labextmcur.col(1) = arma::max(cubelab, 1);
					//Place the cube in the list
					mat clijk(3, 2);
					clijk.col(0) = arma::floor<vec>((labextmcur.col(0) - m_labextm.col(0)) / (m_labrng / (noclpts - 1)));
					clijk.col(1) = arma::floor<vec>((labextmcur.col(1) - m_labextm.col(0)) / (m_labrng / (noclpts - 1)));					
					clijk.elem(arma::find(clijk > noclpts - 2)).fill(noclpts - 2); //Set elements greater than max to max
					std::vector<int> idxlst;
					for (int clidx_k = int(clijk(2,0)); clidx_k <= clijk(2,1); clidx_k++)
					{
						for (int clidx_j = int(clijk(1, 0)); clidx_j <= clijk(1,1); clidx_j++)
						{
							for (int clidx_i = int(clijk(0, 0)); clidx_i <= clijk(0,1); clidx_i++)
							{
								idxlst.push_back(lutInd3D(clidx_i, clidx_j, clidx_k, noclpts-1));
							}
						}
					}
					for (int el = 0; el < idxlst.size(); el++)
					{
						int arr[] = {i, j, k, l};
						try
						{
							m_cubelist[idxlst[el]].push_back(std::vector<int>(arr, arr + 4));
						}
						catch (const std::exception& ex)
						{
							cout << ex.what();
							cout << i << " " << j << " " << k << " " << l << " " << idxlst[el] << "\n";
						}
						//m_cubelist[idxlst[el]].push_back(std::vector<int>(arr, arr + 4));
					}

					//cout << lutInd4D(i, j, k, l, noclpts) << " / " << std::pow(noclpts - 1, 4) << "\n";
				} //i
			} //j
		} //k
		cout << l << "\n";
	} //l
	//Write cubelist to file
	/*std::ofstream file;
	file.open("list.bin", std::ofstream::binary);
	file.write((char*) &nogrdpts, 4);
	file.write((char*) &noclpts, 4);
	for (int i = 0; i < m_cubelist.size(); ++i)
	{
		int size = int(m_cubelist[i].size());
		file.write((char*) &size, 4);
		
		for (int j = 0; j < m_cubelist[i].size(); j++)
		{
			file.write((char*)&m_cubelist[i][j][0], 16);
		}
	}
	file.close();*/
}

RevInt::RevInt(const char* prof, const char* cubelist)
{
	std::ifstream file;
	file.open(cubelist, std::ofstream::binary);
	file.read(reinterpret_cast<char*>(&m_nogrdpts), 4);
	file.read(reinterpret_cast<char*>(&m_noclpts), 4);
	m_cubelist.resize(pow(m_noclpts - 1, 3));
	for (int i=0; i<std::pow(m_noclpts - 1, 3); i++)
	{
		int size;
		file.read(reinterpret_cast<char*>(&size), 4);
		m_cubelist[i].resize(size);
		int arr[4];		
		for (int j=0; j<size; j++)
		{
			file.read(reinterpret_cast<char*>(&arr[0]), 16);
			m_cubelist[i][j] = std::vector<int>(arr, arr + 4);
		}
	}
	file.close();
	m_cfAtoB = makeCForm(prof, "AtoB1");
	m_cfBtoA = makeCForm(prof, "BtoA1");
	m_cmykgrid = permRep(arma::linspace(0, 100, m_nogrdpts), 4, "first");
	m_cmykgrid = m_cmykgrid.t();
	double *labgrid = new double[m_cmykgrid.n_elem / 4 * 3];
	cmsDoTransform(m_cfAtoB, &m_cmykgrid(0, 0), labgrid, m_cmykgrid.n_elem / 4);
	m_labgrid = mat(labgrid, 3, m_cmykgrid.n_elem / 4);
	delete[] labgrid;

	//4D cube vertices
	Mat<int> cv = arma::conv_to<arma::s32_mat>::from(permRep(arma::linspace(0, 1, 2), 4, "first"));

	//4D simplices around main diagonal, 4!=24
	int x[] = { 0, 1, 2, 3 };
	Mat<int> conds = permWoRep(x, 4);
	simplices.resize(24);
	for (int i = 0; i < conds.n_rows; i++)
	{
		Mat<int> set1 = cv.rows(arma::find(cv.col(conds(i, 0)) >= cv.col(conds(i, 1)) == 1));
		Mat<int> set2 = cv.rows(arma::find(cv.col(conds(i, 1)) >= cv.col(conds(i, 2)) == 1));
		Mat<int> set3 = cv.rows(arma::find(cv.col(conds(i, 2)) >= cv.col(conds(i, 3)) == 1));
		Mat<int> set4 = cv.rows(arma::find(cv.col(conds(i, 0)) >= cv.col(conds(i, 2)) == 1));
		Mat<int> set5 = cv.rows(arma::find(cv.col(conds(i, 0)) >= cv.col(conds(i, 3)) == 1));
		Mat<int> set6 = cv.rows(arma::find(cv.col(conds(i, 1)) >= cv.col(conds(i, 3)) == 1));
		Mat<int> set = intersect(intersect(intersect(intersect(intersect(set1, set2), set3), set4), set5), set6);
		simplices[i] = set;
	}

	//Get L*a*b* min and max from grid
	m_labextm.resize(3, 2);
	m_labextm.col(0) = arma::min(m_labgrid, 1);
	m_labextm.col(1) = arma::max(m_labgrid, 1);
	m_labrng = m_labextm.col(1) - m_labextm.col(0);
}

mat RevInt::interp(vec labpt, vec kvals)
{
	//Determine Labpt's incices for cube list and get cubes where it might reside
	vec labptn = labpt - m_labextm.col(0); //Value normalized to minimum
	vec list_ijk = arma::floor(labptn/(m_labrng/(m_noclpts-1)));
	list_ijk(arma::find(list_ijk > m_noclpts - 2)).fill(m_noclpts - 2);
	int clind = lutInd3D(list_ijk(0), list_ijk(1), list_ijk(2), m_noclpts-1);
	
	//Structures to store solutions
	std::vector<arma::s32_mat> sol_simp; //Simplices in which labpt was solved for
	std::vector<std::vector<double> > sol_k;
	std::vector<std::vector<double> > sol_c_rng;
	//cout << clind << " " << m_cubelist[clind].size() << "\n";
	for (int cube = 0; cube < m_cubelist[clind].size(); cube++)
	{		
		std::vector<Mat<int> > simplices_cur = simplices; //Current simplices
		//Search for solutions within current cube's simplices
		for (int simp = 0; simp < 24; simp++)
		{
			simplices_cur[simp].col(0) = simplices_cur[simp].col(0) + m_cubelist[clind][cube][0];
			simplices_cur[simp].col(1) = simplices_cur[simp].col(1) + m_cubelist[clind][cube][1];
			simplices_cur[simp].col(2) = simplices_cur[simp].col(2) + m_cubelist[clind][cube][2];
			simplices_cur[simp].col(3) = simplices_cur[simp].col(3) + m_cubelist[clind][cube][3];

			mat lab_v(3, 5);
			mat cmyk_v(4, 5);
			for (int j = 0; j < 5; j++)
			{
				lab_v.col(j) = m_labgrid.col(lutInd4D(simplices_cur[simp](j, 0), simplices_cur[simp](j, 1), simplices_cur[simp](j, 2), simplices_cur[simp](j, 3), m_nogrdpts));
				cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(simplices_cur[simp](j, 0), simplices_cur[simp](j, 1), simplices_cur[simp](j, 2), simplices_cur[simp](j, 3), m_nogrdpts));
			}

			//SVD solution
			lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
			lab_v.row(lab_v.n_rows - 1) = lab_v.row(0) + lab_v.row(1);

			mat U;
			vec s;
			mat V;
			mat lab_v_svd(lab_v.n_rows, 4);
			lab_v_svd.col(0) = lab_v.col(1) - lab_v.col(0);
			lab_v_svd.col(1) = lab_v.col(2) - lab_v.col(0);
			lab_v_svd.col(2) = lab_v.col(3) - lab_v.col(0);
			lab_v_svd.col(3) = lab_v.col(4) - lab_v.col(0);
			arma::svd(U, s, V, lab_v_svd);
			vec labpt_svd(labpt.n_elem+1);
			labpt_svd(span(0, 2)) = labpt;
			labpt_svd(3) = labpt(0) + labpt(1);
			vec b = labpt_svd - lab_v.col(0);
			s.elem(arma::find(s < 1.0e-13)).zeros();
			vec s_inv = 1 / s;
			s_inv.elem(arma::find_nonfinite(s_inv)).zeros();
			vec bp1 = U.t() * b;
			vec xp = V * diagmat(s_inv) * bp1;
			vec xn = -V.col(3); //Not sure why, SVD returns negative vector (unlike SVD in Matlab)
			vec c0(5);
			vec c1(5);
			c0(span(0, 3)) = -xp / xn;
			c1(span(0, 3)) = (1 - xp) / xn;
			c0(4) = -(sum(xp)) / (sum(xn));
			c1(4) = (1-sum(xp)) / (sum(xn));
			vec c_rng(2);
			c_rng = overlap(c0, c1);

			if (c_rng.is_empty()==false) //If interval width is not zero
			{
				vec xcomplete0 = xp + c_rng(0) * xn;
				vec xcomplete1 = xp + ((c_rng(1)-c_rng(0)) + c_rng(0)) * xn;
				mat simpvecs(cmyk_v.n_rows, 4);
				simpvecs.col(0) = cmyk_v.col(1)-cmyk_v.col(0);
				simpvecs.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				simpvecs.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				simpvecs.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec sol0 = simpvecs * xcomplete0 + cmyk_v.col(0);
				vec sol1 = simpvecs * xcomplete1 + cmyk_v.col(0);
				sol_simp.push_back(simplices_cur[simp]);
				double arr_k[] = {sol0(3), sol1(3)};
				sol_k.push_back(std::vector<double>(arr_k, arr_k+2));
				double arr_c_rng[] = {c_rng(0), c_rng(1)};
				sol_c_rng.push_back(std::vector<double>(arr_c_rng, arr_c_rng + 2));
			}
		} //simp
	} //cube

	//If no solutions are found
	if (sol_simp.size()==0)
	{
		mat res;
		return res;
	}

	//Clip K values to range

	double kmin = sol_k[0][0];
	double kmax = sol_k[0][0];
	for (int i = 0; i < sol_k.size(); i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (sol_k[i][j] < kmin)
			{
				kmin = sol_k[i][j];
			}
			if (sol_k[i][j] > kmax)
			{
				kmax = sol_k[i][j];
			}
		}
	}
	kvals(arma::find(kvals < kmin)).fill(kmin);
	kvals(arma::find(kvals > kmax)).fill(kmax);
	//Solve for desired K values
	mat res(4, kvals.n_elem);
	res.fill(arma::datum::nan);
	//Get unique kvals (black ink amounts)
	uvec unqkvals = arma::find_unique(kvals);
	uvec::iterator it = unqkvals.begin();
	uvec::iterator it_end = unqkvals.end();
	for (; it!=it_end; it++)
	{
		//Get current kval from unique values list
		double kval = kvals(*it);
		//Get indices where current kval occurs
		uvec kvalidx = arma::find(kvals==kval);
		for (int rngidx = 0; rngidx < sol_k.size(); rngidx++)
		{
			if (kval >= std::min(sol_k[rngidx][0], sol_k[rngidx][1]) && kval <= std::max(sol_k[rngidx][0], sol_k[rngidx][1]))
			{
				//kfct factor between min and max K within one simplex
				double kfct = 0;
				if (sol_k[rngidx][0] != sol_k[rngidx][1]) //Avoid division by zero if min=max in current simplex
				{
					kfct = (kval - std::min(sol_k[rngidx][0], sol_k[rngidx][1])) / (std::max(sol_k[rngidx][0], sol_k[rngidx][1]) - std::min(sol_k[rngidx][0], sol_k[rngidx][1]));

				}
				
				mat lab_v(3, 5);
				mat cmyk_v(4, 5);
				for (int j = 0; j < 5; j++)
				{
					lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
					cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
				}

				//SVD solution
				lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
				lab_v.row(lab_v.n_rows - 1) = lab_v.row(0) + lab_v.row(1);

				mat U;
				vec s;
				mat V;
				mat lab_v_svd(lab_v.n_rows, 4);
				lab_v_svd.col(0) = lab_v.col(1) - lab_v.col(0);
				lab_v_svd.col(1) = lab_v.col(2) - lab_v.col(0);
				lab_v_svd.col(2) = lab_v.col(3) - lab_v.col(0);
				lab_v_svd.col(3) = lab_v.col(4) - lab_v.col(0);
				arma::svd(U, s, V, lab_v_svd);
				vec labpt_svd(labpt.n_elem + 1);
				labpt_svd(span(0, 2)) = labpt;
				labpt_svd(3) = labpt(0) + labpt(1);
				vec b = labpt_svd - lab_v.col(0);
				s.elem(arma::find(s < 1.0e-13)).zeros();
				vec s_inv = 1 / s;
				s_inv.elem(arma::find_nonfinite(s_inv)).zeros();
				vec bp1 = U.t() * b;
				vec xp = V * diagmat(s_inv) * bp1;
				vec xn = -V.col(3); //Not sure why, SVD returns negative vector (unlike SVD in Matlab)
				
				//Does K increase or decrease as c changes?
				if (sol_k[rngidx][0] >= sol_k[rngidx][1])
				{
					kfct = 1 - kfct;
				}
				vec xcomplete = xp + ((sol_c_rng[rngidx][1] - sol_c_rng[rngidx][0]) * kfct + sol_c_rng[rngidx][0]) * xn;
				mat cmyksol(4, 4);
				cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec solution = cmyksol * xcomplete + cmyk_v.col(0);
				res.each_col(kvalidx) = solution;
				break;
			}
		}
	}
	return res;
}

mat RevInt::interpex(vec labpt, vec kvals)
{
	//Determine Labpt's incices for cube list and get cubes where it might reside
	vec labptn = labpt - m_labextm.col(0); //Value normalized to minimum
	vec list_ijk = arma::floor(labptn / (m_labrng / (m_noclpts - 1)));
	list_ijk(arma::find(list_ijk > m_noclpts - 2)).fill(m_noclpts - 2);
	int clind = lutInd3D(list_ijk(0), list_ijk(1), list_ijk(2), m_noclpts - 1);

	//Structures to store solutions
	std::vector<arma::s32_mat> sol_simp; //Simplices in which labpt was solved for
	std::vector<std::vector<double> > sol_k;
	std::vector<std::vector<double> > sol_c_rng;
	//cout << clind << " " << m_cubelist[clind].size() << "\n";
	for (int cube = 0; cube < m_cubelist[clind].size(); cube++)
	{
		std::vector<Mat<int> > simplices_cur = simplices; //Current simplices
		//Search for solutions within current cube's simplices
		for (int simp = 0; simp < 24; simp++)
		{
			simplices_cur[simp].col(0) = simplices_cur[simp].col(0) + m_cubelist[clind][cube][0];
			simplices_cur[simp].col(1) = simplices_cur[simp].col(1) + m_cubelist[clind][cube][1];
			simplices_cur[simp].col(2) = simplices_cur[simp].col(2) + m_cubelist[clind][cube][2];
			simplices_cur[simp].col(3) = simplices_cur[simp].col(3) + m_cubelist[clind][cube][3];

			mat lab_v(3, 5);
			mat cmyk_v(4, 5);
			for (int j = 0; j < 5; j++)
			{
				lab_v.col(j) = m_labgrid.col(lutInd4D(simplices_cur[simp](j, 0), simplices_cur[simp](j, 1), simplices_cur[simp](j, 2), simplices_cur[simp](j, 3), m_nogrdpts));
				cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(simplices_cur[simp](j, 0), simplices_cur[simp](j, 1), simplices_cur[simp](j, 2), simplices_cur[simp](j, 3), m_nogrdpts));
			}

			//SVD solution
			lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
			lab_v.row(lab_v.n_rows - 1) = lab_v.row(0) + lab_v.row(1);

			mat U;
			vec s;
			mat V;
			mat lab_v_svd(lab_v.n_rows, 4);
			lab_v_svd.col(0) = lab_v.col(1) - lab_v.col(0);
			lab_v_svd.col(1) = lab_v.col(2) - lab_v.col(0);
			lab_v_svd.col(2) = lab_v.col(3) - lab_v.col(0);
			lab_v_svd.col(3) = lab_v.col(4) - lab_v.col(0);
			arma::svd(U, s, V, lab_v_svd);
			vec labpt_svd(labpt.n_elem + 1);
			labpt_svd(span(0, 2)) = labpt;
			labpt_svd(3) = labpt(0) + labpt(1);
			vec b = labpt_svd - lab_v.col(0);
			s.elem(arma::find(s < 1.0e-13)).zeros();
			vec s_inv = 1 / s;
			s_inv.elem(arma::find_nonfinite(s_inv)).zeros();
			vec bp1 = U.t() * b;
			vec xp = V * diagmat(s_inv) * bp1;
			vec xn = -V.col(3); //Not sure why, SVD returns negative vector (unlike SVD in Matlab)
			vec c0(5);
			vec c1(5);
			c0(span(0, 3)) = -xp / xn;
			c1(span(0, 3)) = (1 - xp) / xn;
			c0(4) = -(sum(xp)) / (sum(xn));
			c1(4) = (1 - sum(xp)) / (sum(xn));
			vec c_rng(2);
			c_rng = overlap(c0, c1);

			if (c_rng.is_empty() == false) //If interval width is not zero
			{
				vec xcomplete0 = xp + c_rng(0) * xn;
				vec xcomplete1 = xp + ((c_rng(1) - c_rng(0)) + c_rng(0)) * xn;
				mat simpvecs(cmyk_v.n_rows, 4);
				simpvecs.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				simpvecs.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				simpvecs.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				simpvecs.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec sol0 = simpvecs * xcomplete0 + cmyk_v.col(0);
				vec sol1 = simpvecs * xcomplete1 + cmyk_v.col(0);
				sol_simp.push_back(simplices_cur[simp]);
				double arr_k[] = { sol0(3), sol1(3) };
				sol_k.push_back(std::vector<double>(arr_k, arr_k + 2));
				double arr_c_rng[] = { c_rng(0), c_rng(1) };
				sol_c_rng.push_back(std::vector<double>(arr_c_rng, arr_c_rng + 2));
			}
		} //simp
	} //cube

	//If no solutions are found
	if (sol_simp.size() == 0)
	{
		mat res;
		return res;
	}

	//Clip K values to range

	double kmin = sol_k[0][0];
	double kmax = sol_k[0][0];
	for (int i = 0; i < sol_k.size(); i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (sol_k[i][j] < kmin)
			{
				kmin = sol_k[i][j];
			}
			if (sol_k[i][j] > kmax)
			{
				kmax = sol_k[i][j];
			}
		}
	}
	uvec klows = arma::find(kvals < kmin);
	uvec khis = arma::find(kvals > kmax);
	if (klows.n_elem > 0)
	{
		kvals(klows).fill(kvals(klows(klows.n_elem - 1)));
	}
	if (khis.n_elem > 0)
	{
		kvals(khis).fill(kvals(khis(0)));
	}

	//Solve for desired K values
	mat res(4, kvals.n_elem);
	res.fill(arma::datum::nan);

	//If some values in kvals are less than kmin
	if (klows.n_elem > 0)
	{
		//Find simplex which solves for kmin and use it to extrapolate for K's below kmin
		for (int rngidx = 0; rngidx < sol_k.size(); rngidx++)
		{
			if (kmin >= std::min(sol_k[rngidx][0], sol_k[rngidx][1]) && kmin <= std::max(sol_k[rngidx][0], sol_k[rngidx][1]))
			{
				mat lab_v(3, 5);
				mat cmyk_v(4, 5);
				for (int j = 0; j < 5; j++)
				{
					lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
					cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
				}

				lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
				lab_v.row(lab_v.n_rows - 1) = cmyk_v.row(cmyk_v.n_rows - 1);
								
				mat lab_v_sol(lab_v.n_rows, 4);
				lab_v_sol.col(0) = lab_v.col(1) - lab_v.col(0);
				lab_v_sol.col(1) = lab_v.col(2) - lab_v.col(0);
				lab_v_sol.col(2) = lab_v.col(3) - lab_v.col(0);
				lab_v_sol.col(3) = lab_v.col(4) - lab_v.col(0);
				vec labptk(4);
				labptk(span(0, 2)) = labpt;
				labptk(3) = kvals(0);
				vec b = labptk - lab_v.col(0);
				vec sol = arma::solve(lab_v_sol, b);
				mat cmyksol(4, 4);
				cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec solution = cmyksol * sol + cmyk_v.col(0);
				res.each_col(klows) = solution;

				break;						
			}
		}
	}

	//If some values in kvals are greater than kmax
	if (khis.n_elem > 0)
	{
		//Find simplex which solves for kmax and use it to extrapolate for K's above kmin
		for (int rngidx = 0; rngidx < sol_k.size(); rngidx++)
		{
			if (kmax >= std::min(sol_k[rngidx][0], sol_k[rngidx][1]) && kmax <= std::max(sol_k[rngidx][0], sol_k[rngidx][1]))
			{
				mat lab_v(3, 5);
				mat cmyk_v(4, 5);
				for (int j = 0; j < 5; j++)
				{
					lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
					cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
				}

				lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
				lab_v.row(lab_v.n_rows - 1) = cmyk_v.row(cmyk_v.n_rows - 1);

				mat lab_v_sol(lab_v.n_rows, 4);
				lab_v_sol.col(0) = lab_v.col(1) - lab_v.col(0);
				lab_v_sol.col(1) = lab_v.col(2) - lab_v.col(0);
				lab_v_sol.col(2) = lab_v.col(3) - lab_v.col(0);
				lab_v_sol.col(3) = lab_v.col(4) - lab_v.col(0);
				vec labptk(4);
				labptk(span(0, 2)) = labpt;
				labptk(3) = kvals(kvals.n_elem - 1);
				vec b = labptk - lab_v.col(0);
				vec sol = arma::solve(lab_v_sol, b);
				mat cmyksol(4, 4);
				cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec solution = cmyksol * sol + cmyk_v.col(0);
				res.each_col(khis) = solution;

				break;
			}
		}
	}
	
	//Solve for K's within [kmin, kmax]
	vec::iterator it;
	//Set iterator lower bound (K's in kvals greater than kmin)
	if (klows.n_elem > 0)
	{
		it = &kvals(klows(klows.n_elem - 1));
		it++;
	}
	else
	{
		it = kvals.begin();
	}
	//Set iterator upper bound (K's in kvals less than kmax)
	vec::iterator it_end = nullptr;
	if (khis.n_elem > 0)
	{
		it_end = &kvals(khis(0));
	}
	else
	{
		it_end = kvals.end();
	}

	for (; it != it_end; it++)
	{
		//Get current kval from unique values list
		double kval = *it;
		//Get indices where current kval occurs
		uvec kvalidx = arma::find(kvals == kval);
		for (int rngidx = 0; rngidx < sol_k.size(); rngidx++)
		{
			if (kval >= std::min(sol_k[rngidx][0], sol_k[rngidx][1]) && kval <= std::max(sol_k[rngidx][0], sol_k[rngidx][1]))
			{
				//kfct factor between min and max K within one simplex
				double kfct = 0;
				if (sol_k[rngidx][0] != sol_k[rngidx][1]) //Avoid division by zero if min=max in current simplex
				{
					kfct = (kval - std::min(sol_k[rngidx][0], sol_k[rngidx][1])) / (std::max(sol_k[rngidx][0], sol_k[rngidx][1]) - std::min(sol_k[rngidx][0], sol_k[rngidx][1]));

				}

				mat lab_v(3, 5);
				mat cmyk_v(4, 5);
				for (int j = 0; j < 5; j++)
				{
					lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
					cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
				}

				//SVD solution
				lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
				lab_v.row(lab_v.n_rows - 1) = lab_v.row(0) + lab_v.row(1);

				mat U;
				vec s;
				mat V;
				mat lab_v_svd(lab_v.n_rows, 4);
				lab_v_svd.col(0) = lab_v.col(1) - lab_v.col(0);
				lab_v_svd.col(1) = lab_v.col(2) - lab_v.col(0);
				lab_v_svd.col(2) = lab_v.col(3) - lab_v.col(0);
				lab_v_svd.col(3) = lab_v.col(4) - lab_v.col(0);
				arma::svd(U, s, V, lab_v_svd);
				vec labpt_svd(labpt.n_elem + 1);
				labpt_svd(span(0, 2)) = labpt;
				labpt_svd(3) = labpt(0) + labpt(1);
				vec b = labpt_svd - lab_v.col(0);
				s.elem(arma::find(s < 1.0e-13)).zeros();
				vec s_inv = 1 / s;
				s_inv.elem(arma::find_nonfinite(s_inv)).zeros();
				vec bp1 = U.t() * b;
				vec xp = V * diagmat(s_inv) * bp1;
				vec xn = -V.col(3); //Not sure why, SVD returns negative vector (unlike SVD in Matlab)

				//Does K increase or decrease as c changes?
				if (sol_k[rngidx][0] >= sol_k[rngidx][1])
				{
					kfct = 1 - kfct;
				}
				vec xcomplete = xp + ((sol_c_rng[rngidx][1] - sol_c_rng[rngidx][0]) * kfct + sol_c_rng[rngidx][0]) * xn;
				mat cmyksol(4, 4);
				cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec solution = cmyksol * xcomplete + cmyk_v.col(0);
				res.each_col(kvalidx) = solution;
				break;
			}
		}
	}
	return res;
}

mat RevInt::interpexLut(vec labpt, vec kvals, std::vector<arma::s32_mat> *sol_simp_vec,
	std::vector<std::vector<double> > *sol_k_vec, double *kminout, double *kmaxout)
{
	//Determine Labpt's incices for cube list and get cubes where it might reside
	vec labptn = labpt - m_labextm.col(0); //Value normalized to minimum
	vec list_ijk = arma::floor(labptn / (m_labrng / (m_noclpts - 1)));
	list_ijk(arma::find(list_ijk > m_noclpts - 2)).fill(m_noclpts - 2);
	int clind = lutInd3D(list_ijk(0), list_ijk(1), list_ijk(2), m_noclpts - 1);

	//Structures to store solutions
	std::vector<arma::s32_mat> sol_simp; //Simplices in which labpt was solved for
	std::vector<std::vector<double> > sol_k;
	std::vector<std::vector<double> > sol_c_rng;
	for (int cube = 0; cube < m_cubelist[clind].size(); cube++)
	{
		std::vector<Mat<int> > simplices_cur = simplices; //Current simplices
		//Search for solutions within current cube's simplices
		for (int simp = 0; simp < 24; simp++)
		{
			simplices_cur[simp].col(0) = simplices_cur[simp].col(0) + m_cubelist[clind][cube][0];
			simplices_cur[simp].col(1) = simplices_cur[simp].col(1) + m_cubelist[clind][cube][1];
			simplices_cur[simp].col(2) = simplices_cur[simp].col(2) + m_cubelist[clind][cube][2];
			simplices_cur[simp].col(3) = simplices_cur[simp].col(3) + m_cubelist[clind][cube][3];

			mat lab_v(3, 5);
			mat cmyk_v(4, 5);
			for (int j = 0; j < 5; j++)
			{
				lab_v.col(j) = m_labgrid.col(lutInd4D(simplices_cur[simp](j, 0), simplices_cur[simp](j, 1), simplices_cur[simp](j, 2), simplices_cur[simp](j, 3), m_nogrdpts));
				cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(simplices_cur[simp](j, 0), simplices_cur[simp](j, 1), simplices_cur[simp](j, 2), simplices_cur[simp](j, 3), m_nogrdpts));
			}

			//SVD solution
			lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
			lab_v.row(lab_v.n_rows - 1) = lab_v.row(0) + lab_v.row(1);

			mat U;
			vec s;
			mat V;
			mat lab_v_svd(lab_v.n_rows, 4);
			lab_v_svd.col(0) = lab_v.col(1) - lab_v.col(0);
			lab_v_svd.col(1) = lab_v.col(2) - lab_v.col(0);
			lab_v_svd.col(2) = lab_v.col(3) - lab_v.col(0);
			lab_v_svd.col(3) = lab_v.col(4) - lab_v.col(0);
			arma::svd(U, s, V, lab_v_svd);
			vec labpt_svd(labpt.n_elem + 1);
			labpt_svd(span(0, 2)) = labpt;
			labpt_svd(3) = labpt(0) + labpt(1);
			vec b = labpt_svd - lab_v.col(0);
			s.elem(arma::find(s < 1.0e-13)).zeros();
			vec s_inv = 1 / s;
			s_inv.elem(arma::find_nonfinite(s_inv)).zeros();
			vec bp1 = U.t() * b;
			vec xp = V * diagmat(s_inv) * bp1;
			vec xn = -V.col(3); //Not sure why, SVD returns negative vector (unlike SVD in Matlab)
			vec c0(5);
			vec c1(5);
			c0(span(0, 3)) = -xp / xn;
			c1(span(0, 3)) = (1 - xp) / xn;
			c0(4) = -(sum(xp)) / (sum(xn));
			c1(4) = (1 - sum(xp)) / (sum(xn));
			vec c_rng(2);
			c_rng = overlap(c0, c1);

			if (c_rng.is_empty() == false) //If interval width is not zero
			{
				vec xcomplete0 = xp + c_rng(0) * xn;
				vec xcomplete1 = xp + ((c_rng(1) - c_rng(0)) + c_rng(0)) * xn;
				mat simpvecs(cmyk_v.n_rows, 4);
				simpvecs.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				simpvecs.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				simpvecs.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				simpvecs.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec sol0 = simpvecs * xcomplete0 + cmyk_v.col(0);
				vec sol1 = simpvecs * xcomplete1 + cmyk_v.col(0);
				sol_simp.push_back(simplices_cur[simp]);
				double arr_k[] = { sol0(3), sol1(3) };
				sol_k.push_back(std::vector<double>(arr_k, arr_k + 2));
				double arr_c_rng[] = { c_rng(0), c_rng(1) };
				sol_c_rng.push_back(std::vector<double>(arr_c_rng, arr_c_rng + 2));
			}
		} //simp
	} //cube

	//If no solutions are found
	if (sol_simp.size() == 0)
	{
		mat res;
		return res;
	}
	

	//Clip K values to range

	double kmin = sol_k[0][0];
	double kmax = sol_k[0][0];
	for (int i = 0; i < sol_k.size(); i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (sol_k[i][j] < kmin)
			{
				kmin = sol_k[i][j];
			}
			if (sol_k[i][j] > kmax)
			{
				kmax = sol_k[i][j];
			}
		}
	}

	//Update structures in caller function
	*kminout = kmin;
	*kmaxout = kmax;
	*sol_simp_vec = sol_simp;
	*sol_k_vec = sol_k;

	uvec klows = arma::find(kvals < kmin);
	uvec khis = arma::find(kvals > kmax);
	if (klows.n_elem > 0)
	{
		kvals(klows).fill(kvals(klows(klows.n_elem - 1)));
	}
	if (khis.n_elem > 0)
	{
		kvals(khis).fill(kvals(khis(0)));
	}

	//Solve for desired K values
	mat res(4, kvals.n_elem);
	res.fill(arma::datum::nan);

	//If some values in kvals are less than kmin
	if (klows.n_elem > 0)
	{
		//Find simplex which solves for kmin and use it to extrapolate for K's below kmin
		for (int rngidx = 0; rngidx < sol_k.size(); rngidx++)
		{
			if (kmin >= std::min(sol_k[rngidx][0], sol_k[rngidx][1]) && kmin <= std::max(sol_k[rngidx][0], sol_k[rngidx][1]))
			{
				mat lab_v(3, 5);
				mat cmyk_v(4, 5);
				for (int j = 0; j < 5; j++)
				{
					lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
					cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
				}

				lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
				lab_v.row(lab_v.n_rows - 1) = cmyk_v.row(cmyk_v.n_rows - 1);

				mat lab_v_sol(lab_v.n_rows, 4);
				lab_v_sol.col(0) = lab_v.col(1) - lab_v.col(0);
				lab_v_sol.col(1) = lab_v.col(2) - lab_v.col(0);
				lab_v_sol.col(2) = lab_v.col(3) - lab_v.col(0);
				lab_v_sol.col(3) = lab_v.col(4) - lab_v.col(0);
				vec labptk(4);
				labptk(span(0, 2)) = labpt;
				labptk(3) = kvals(0);
				vec b = labptk - lab_v.col(0);
				vec sol = arma::solve(lab_v_sol, b);
				mat cmyksol(4, 4);
				cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec solution = cmyksol * sol + cmyk_v.col(0);
				res.each_col(klows) = solution;

				break;
			}
		}
	}

	//If some values in kvals are greater than kmax
	if (khis.n_elem > 0)
	{
		//Find simplex which solves for kmax and use it to extrapolate for K's above kmin
		for (int rngidx = 0; rngidx < sol_k.size(); rngidx++)
		{
			if (kmax >= std::min(sol_k[rngidx][0], sol_k[rngidx][1]) && kmax <= std::max(sol_k[rngidx][0], sol_k[rngidx][1]))
			{
				mat lab_v(3, 5);
				mat cmyk_v(4, 5);
				for (int j = 0; j < 5; j++)
				{
					lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
					cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
				}

				lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
				lab_v.row(lab_v.n_rows - 1) = cmyk_v.row(cmyk_v.n_rows - 1);

				mat lab_v_sol(lab_v.n_rows, 4);
				lab_v_sol.col(0) = lab_v.col(1) - lab_v.col(0);
				lab_v_sol.col(1) = lab_v.col(2) - lab_v.col(0);
				lab_v_sol.col(2) = lab_v.col(3) - lab_v.col(0);
				lab_v_sol.col(3) = lab_v.col(4) - lab_v.col(0);
				vec labptk(4);
				labptk(span(0, 2)) = labpt;
				labptk(3) = kvals(kvals.n_elem - 1);
				vec b = labptk - lab_v.col(0);
				vec sol = arma::solve(lab_v_sol, b);
				mat cmyksol(4, 4);
				cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec solution = cmyksol * sol + cmyk_v.col(0);
				res.each_col(khis) = solution;

				break;
			}
		}
	}

	//Solve for K's within [kmin, kmax]
	vec::iterator it;
	//Set iterator lower bound (K's in kvals greater than kmin)
	if (klows.n_elem > 0)
	{
		it = &kvals(klows(klows.n_elem - 1));
		it++;
	}
	else
	{
		it = kvals.begin();
	}
	//Set iterator upper bound (K's in kvals less than kmax)
	vec::iterator it_end = nullptr;
	if (khis.n_elem > 0)
	{
		it_end = &kvals(khis(0));
	}
	else
	{
		it_end = kvals.end();
	}

	for (; it != it_end; it++)
	{
		//Get current kval from unique values list
		double kval = *it;
		//Get indices where current kval occurs
		uvec kvalidx = arma::find(kvals == kval);
		for (int rngidx = 0; rngidx < sol_k.size(); rngidx++)
		{
			if (kval >= std::min(sol_k[rngidx][0], sol_k[rngidx][1]) && kval <= std::max(sol_k[rngidx][0], sol_k[rngidx][1]))
			{
				//kfct factor between min and max K within one simplex
				double kfct = 0;
				if (sol_k[rngidx][0] != sol_k[rngidx][1]) //Avoid division by zero if min=max in current simplex
				{
					kfct = (kval - std::min(sol_k[rngidx][0], sol_k[rngidx][1])) / (std::max(sol_k[rngidx][0], sol_k[rngidx][1]) - std::min(sol_k[rngidx][0], sol_k[rngidx][1]));

				}

				mat lab_v(3, 5);
				mat cmyk_v(4, 5);
				for (int j = 0; j < 5; j++)
				{
					lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
					cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp[rngidx](j, 0), sol_simp[rngidx](j, 1), sol_simp[rngidx](j, 2), sol_simp[rngidx](j, 3), m_nogrdpts));
				}

				//SVD solution
				lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
				lab_v.row(lab_v.n_rows - 1) = lab_v.row(0) + lab_v.row(1);

				mat U;
				vec s;
				mat V;
				mat lab_v_svd(lab_v.n_rows, 4);
				lab_v_svd.col(0) = lab_v.col(1) - lab_v.col(0);
				lab_v_svd.col(1) = lab_v.col(2) - lab_v.col(0);
				lab_v_svd.col(2) = lab_v.col(3) - lab_v.col(0);
				lab_v_svd.col(3) = lab_v.col(4) - lab_v.col(0);
				arma::svd(U, s, V, lab_v_svd);
				vec labpt_svd(labpt.n_elem + 1);
				labpt_svd(span(0, 2)) = labpt;
				labpt_svd(3) = labpt(0) + labpt(1);
				vec b = labpt_svd - lab_v.col(0);
				s.elem(arma::find(s < 1.0e-13)).zeros();
				vec s_inv = 1 / s;
				s_inv.elem(arma::find_nonfinite(s_inv)).zeros();
				vec bp1 = U.t() * b;
				vec xp = V * diagmat(s_inv) * bp1;
				vec xn = -V.col(3); //Not sure why, SVD returns negative vector (unlike SVD in Matlab)

				//Does K increase or decrease as c changes?
				if (sol_k[rngidx][0] >= sol_k[rngidx][1])
				{
					kfct = 1 - kfct;
				}
				vec xcomplete = xp + ((sol_c_rng[rngidx][1] - sol_c_rng[rngidx][0]) * kfct + sol_c_rng[rngidx][0]) * xn;
				mat cmyksol(4, 4);
				cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
				cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
				cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
				cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
				vec solution = cmyksol * xcomplete + cmyk_v.col(0);
				res.each_col(kvalidx) = solution;
				break;
			}
		}
	}
	return res;
}

void RevInt::interpLut(double *lutcmykarr, unsigned short *outgamarr, double* kminmax, int nogrdpts)
{
	mat lutlab = permRep(arma::linspace(0, 1, nogrdpts), 4, "last");
	lutlab = lutlab.t();
	lutlab.each_col() %= m_labrng;
	lutlab.each_col() += m_labextm.col(0);
	mat lutcmyk(lutcmykarr, lutlab.n_rows, lutlab.n_cols, false, true);
	lutcmyk.fill(0);
	arma::Col<arma::u16> gmaoutgam(outgamarr, std::pow(nogrdpts, 3));
	arma::Col<arma::u16> solved(gmaoutgam.n_elem);
	solved.fill(0);
	mat kminmaxmat(kminmax, 2, gmaoutgam.n_elem, false, true);
	kminmaxmat.fill(0);

	for (int l = 0; l < nogrdpts - 1; l++)
	{
		for (int k = 0; k < nogrdpts - 1; k++)
		{
			for (int j = 0; j < nogrdpts - 1; j++)
			{
				//Which Lab points are in-gamut and which are out-of-gamut?
				arma::Col<arma::u16> oog(8);
				uvec inds3d(8);
				uvec inds4d(8);
				inds3d(0) = lutInd3D(l, k, j, nogrdpts);
				inds4d(0) = lutInd4D(l, k, j, 0, nogrdpts);
				oog(0) = gmaoutgam(inds3d(0));

				inds3d(1) = lutInd3D(l, k, j+1, nogrdpts);
				inds4d(1) = lutInd4D(l, k, j+1, 0, nogrdpts);
				oog(1) = gmaoutgam(inds3d(1));

				inds3d(2) = lutInd3D(l, k+1, j, nogrdpts);
				inds4d(2) = lutInd4D(l, k+1, j, 0, nogrdpts);
				oog(2) = gmaoutgam(inds3d(2));

				inds3d(3) = lutInd3D(l, k+1, j+1, nogrdpts);
				inds4d(3) = lutInd4D(l, k+1, j+1, 0, nogrdpts);
				oog(3) = gmaoutgam(inds3d(3));

				inds3d(4) = lutInd3D(l+1, k, j, nogrdpts);
				inds4d(4) = lutInd4D(l+1, k, j, 0, nogrdpts);
				oog(4) = gmaoutgam(inds3d(4));

				inds3d(5) = lutInd3D(l+1, k, j+1, nogrdpts);
				inds4d(5) = lutInd4D(l+1, k, j+1, 0, nogrdpts);
				oog(5) = gmaoutgam(inds3d(5));

				inds3d(6) = lutInd3D(l+1, k+1, j, nogrdpts);
				inds4d(6) = lutInd4D(l+1, k+1, j, 0, nogrdpts);
				oog(6) = gmaoutgam(inds3d(6));

				inds3d(7) = lutInd3D(l+1, k+1, j+1, nogrdpts);
				inds4d(7) = lutInd4D(l+1, k+1, j+1, 0, nogrdpts);
				oog(7) = gmaoutgam(inds3d(7));

				//If not all Lab points are out of gamut
				if (arma::sum(oog) != 8)
				{
					//Structures to hold solution for in-gamut points
					std::vector<std::vector<arma::s32_mat> > sol_simp_vec(8); //Simplices in which labpt was solved for
					std::vector<std::vector<std::vector<double> > > sol_k_vec(8);
					double kmins[8];
					double kmaxs[8];

					//Solve for in-gamut Lab points first
					uvec ingam = arma::find(oog == 0);
					for (uvec::iterator it = ingam.begin(); it != ingam.end(); it++)
					{
						//Continue if point is already solved
						if (solved(inds3d(*it)) == 1)
						{
							continue;
						}
						//Call method which returns cmyk solutions to labpt and updates solution structures and Kmin Kmax
						mat A = interpexLut(lutlab.col(inds4d(*it)), arma::linspace(0, 100, nogrdpts),
							&sol_simp_vec[*it], &sol_k_vec[*it], &kmins[*it], &kmaxs[*it]);
						if (A.is_empty()) //If there is no solution although point appears to be in-gamut
						{
							oog(*it) = 1;
							gmaoutgam(inds3d(*it)) = 1;
							continue;
						}
						if (A.has_nan())
						{
							cout << *it << " " << sol_k_vec[*it].size() << "\n";
							oog(*it) = 1;
							gmaoutgam(inds3d(*it)) = 1;
							continue;
						}

						//Previous version, goes wrong if no solution is found
						/*lutcmyk(span(0, 3), span(inds4d(*it), inds4d(*it) + nogrdpts - 1)) =
							interpexLut(lutlab.col(inds4d(*it)), arma::linspace(0, 100, nogrdpts),
								&sol_simp_vec[*it], &sol_k_vec[*it], &kmins[*it], &kmaxs[*it]);*/
						
						//Place in main matrix
						lutcmyk(span(0, 3), span(inds4d(*it), inds4d(*it) + nogrdpts - 1)) = A;
						//A.t().print();
						//Mark vertex as solved
						solved(inds3d(*it)) = 1;
						//Fill kmin kmax matrix
						kminmaxmat(0, inds3d(*it)) = kmins[*it];
						kminmaxmat(1, inds3d(*it)) = kmaxs[*it];
					}

					//Solve for out-of-gamut Lab points
					uvec outgam = arma::find(oog == 1);
					if (outgam.n_elem == 8) //If all points are out-of gamut although some appeared to be in-gamut
					{
						continue;
					}
					ingam = arma::find(oog == 0); //Update list of in-gamut points

					//If all in-gamut points were solved in previous iterations, solve for one
					//to update solution strucures to be used for out-of gamut points
					size_t sizesum = 0;
					for (int solk = 0; solk < 8; solk++)
					{
						sizesum += sol_k_vec[solk].size();
					}
					//if (sizesum == 0)
					//{
						//Call method which returns cmyk solutions to labpt and updates solution structures and Kmin Kmax
						mat A = interpexLut(lutlab.col(inds4d(ingam(0))), arma::linspace(0, 100, nogrdpts),
							&sol_simp_vec[ingam(0)], &sol_k_vec[ingam(0)], &kmins[ingam(0)], &kmaxs[ingam(0)]);
					//}
					
					for (uvec::iterator it = outgam.begin(); it != outgam.end(); it++)
					{
						//Continue if point is already solved
						if (solved(inds3d(*it)) == 1)
						{
							continue;
						}
						int refptlstind = ingam(0); //Index to reference point (any in-gamut) in list of 8 points
						//cout << sol_k_vec[refptlstind].size() << "\n";
						vec labpt = lutlab.col(inds4d(*it));
						vec kvals = arma::linspace(0, 100, nogrdpts);
						uvec klows = arma::find(kvals < kmins[refptlstind]);
						uvec khis = arma::find(kvals > kmaxs[refptlstind]);
						if (klows.n_elem > 0)
						{
							kvals(klows).fill(kvals(klows(klows.n_elem - 1)));
						}
						if (khis.n_elem > 0)
						{
							kvals(khis).fill(kvals(khis(0)));
						}

						//Solve for desired K values
						mat res(4, kvals.n_elem);
						res.fill(arma::datum::nan);


						//If some values in kvals are less than kmin
						if (klows.n_elem > 0)
						{
							//Find simplex which solves for kmin and use it to extrapolate for K's below kmin
							for (int rngidx = 0; rngidx < sol_k_vec[refptlstind].size(); rngidx++)
							{
								if (kmins[refptlstind] >= std::min(sol_k_vec[refptlstind][rngidx][0], sol_k_vec[refptlstind][rngidx][1]) && kmins[refptlstind] <= std::max(sol_k_vec[refptlstind][rngidx][0], sol_k_vec[refptlstind][rngidx][1]))
								{
									mat lab_v(3, 5);
									mat cmyk_v(4, 5);
									for (int j = 0; j < 5; j++)
									{
										lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp_vec[refptlstind][rngidx](j, 0), sol_simp_vec[refptlstind][rngidx](j, 1), sol_simp_vec[refptlstind][rngidx](j, 2), sol_simp_vec[refptlstind][rngidx](j, 3), m_nogrdpts));
										cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp_vec[refptlstind][rngidx](j, 0), sol_simp_vec[refptlstind][rngidx](j, 1), sol_simp_vec[refptlstind][rngidx](j, 2), sol_simp_vec[refptlstind][rngidx](j, 3), m_nogrdpts));
									}

									lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
									lab_v.row(lab_v.n_rows - 1) = cmyk_v.row(cmyk_v.n_rows - 1);

									mat lab_v_sol(lab_v.n_rows, 4);
									lab_v_sol.col(0) = lab_v.col(1) - lab_v.col(0);
									lab_v_sol.col(1) = lab_v.col(2) - lab_v.col(0);
									lab_v_sol.col(2) = lab_v.col(3) - lab_v.col(0);
									lab_v_sol.col(3) = lab_v.col(4) - lab_v.col(0);
									vec labptk(4);
									labptk(span(0, 2)) = labpt;
									labptk(3) = kvals(0);
									vec b = labptk - lab_v.col(0);
									vec sol = arma::solve(lab_v_sol, b);
									mat cmyksol(4, 4);
									cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
									cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
									cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
									cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
									vec solution = cmyksol * sol + cmyk_v.col(0);
									res.each_col(klows) = solution;

									if (solution.has_nan())
									{
										cout << "NaN here" << "\n";
									}

									break;
								}
							}
						}

						//If some values in kvals are greater than kmax
						if (khis.n_elem > 0)
						{
							//Find simplex which solves for kmax and use it to extrapolate for K's above kmin
							for (int rngidx = 0; rngidx < sol_k_vec[refptlstind].size(); rngidx++)
							{
								if (kmaxs[refptlstind] >= std::min(sol_k_vec[refptlstind][rngidx][0], sol_k_vec[refptlstind][rngidx][1]) && kmaxs[refptlstind] <= std::max(sol_k_vec[refptlstind][rngidx][0], sol_k_vec[refptlstind][rngidx][1]))
								{
									mat lab_v(3, 5);
									mat cmyk_v(4, 5);
									for (int j = 0; j < 5; j++)
									{
										lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp_vec[refptlstind][rngidx](j, 0), sol_simp_vec[refptlstind][rngidx](j, 1), sol_simp_vec[refptlstind][rngidx](j, 2), sol_simp_vec[refptlstind][rngidx](j, 3), m_nogrdpts));
										cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp_vec[refptlstind][rngidx](j, 0), sol_simp_vec[refptlstind][rngidx](j, 1), sol_simp_vec[refptlstind][rngidx](j, 2), sol_simp_vec[refptlstind][rngidx](j, 3), m_nogrdpts));
									}

									lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
									lab_v.row(lab_v.n_rows - 1) = cmyk_v.row(cmyk_v.n_rows - 1);

									mat lab_v_sol(lab_v.n_rows, 4);
									lab_v_sol.col(0) = lab_v.col(1) - lab_v.col(0);
									lab_v_sol.col(1) = lab_v.col(2) - lab_v.col(0);
									lab_v_sol.col(2) = lab_v.col(3) - lab_v.col(0);
									lab_v_sol.col(3) = lab_v.col(4) - lab_v.col(0);
									vec labptk(4);
									labptk(span(0, 2)) = labpt;
									labptk(3) = kvals(kvals.n_elem - 1);
									vec b = labptk - lab_v.col(0);
									vec sol = arma::solve(lab_v_sol, b);
									mat cmyksol(4, 4);
									cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
									cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
									cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
									cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
									vec solution = cmyksol * sol + cmyk_v.col(0);
									res.each_col(khis) = solution;

									if (solution.has_nan())
									{
										cout << "NaN here" << "\n";
									}

									break;
								}
							}
						}

						//Solve for K's within [kmin, kmax]
						vec::iterator itk;
						//Set iterator lower bound (K's in kvals greater than kmin)
						if (klows.n_elem > 0)
						{
							itk = &kvals(klows(klows.n_elem - 1));
							itk++;
						}
						else
						{
							itk = kvals.begin();
						}
						//Set iterator upper bound (K's in kvals less than kmax)
						vec::iterator itk_end = nullptr;
						if (khis.n_elem > 0)
						{
							itk_end = &kvals(khis(0));
						}
						else
						{
							itk_end = kvals.end();
						}

						for (; itk != itk_end; itk++)
						{
							//Get current kval from unique values list
							double kval = *itk;
							//Get indices where current kval occurs
							uvec kvalidx = arma::find(kvals == kval);
							for (int rngidx = 0; rngidx < sol_k_vec[refptlstind].size(); rngidx++)
							{
								if (kval >= std::min(sol_k_vec[refptlstind][rngidx][0], sol_k_vec[refptlstind][rngidx][1]) && kval <= std::max(sol_k_vec[refptlstind][rngidx][0], sol_k_vec[refptlstind][rngidx][1]))
								{
									mat lab_v(3, 5);
									mat cmyk_v(4, 5);
									for (int j = 0; j < 5; j++)
									{
										lab_v.col(j) = m_labgrid.col(lutInd4D(sol_simp_vec[refptlstind][rngidx](j, 0), sol_simp_vec[refptlstind][rngidx](j, 1), sol_simp_vec[refptlstind][rngidx](j, 2), sol_simp_vec[refptlstind][rngidx](j, 3), m_nogrdpts));
										cmyk_v.col(j) = m_cmykgrid.col(lutInd4D(sol_simp_vec[refptlstind][rngidx](j, 0), sol_simp_vec[refptlstind][rngidx](j, 1), sol_simp_vec[refptlstind][rngidx](j, 2), sol_simp_vec[refptlstind][rngidx](j, 3), m_nogrdpts));
									}

									lab_v.resize(lab_v.n_rows + 1, lab_v.n_cols);
									lab_v.row(lab_v.n_rows - 1) = cmyk_v.row(cmyk_v.n_rows - 1);

									mat lab_v_sol(lab_v.n_rows, 4);
									lab_v_sol.col(0) = lab_v.col(1) - lab_v.col(0);
									lab_v_sol.col(1) = lab_v.col(2) - lab_v.col(0);
									lab_v_sol.col(2) = lab_v.col(3) - lab_v.col(0);
									lab_v_sol.col(3) = lab_v.col(4) - lab_v.col(0);
									vec labptk(4);
									labptk(span(0, 2)) = labpt;
									labptk(3) = kvals(kvals.n_elem - 1);
									vec b = labptk - lab_v.col(0);
									vec sol = arma::solve(lab_v_sol, b);
									mat cmyksol(4, 4);
									cmyksol.col(0) = cmyk_v.col(1) - cmyk_v.col(0);
									cmyksol.col(1) = cmyk_v.col(2) - cmyk_v.col(0);
									cmyksol.col(2) = cmyk_v.col(3) - cmyk_v.col(0);
									cmyksol.col(3) = cmyk_v.col(4) - cmyk_v.col(0);
									vec solution = cmyksol * sol + cmyk_v.col(0);
									res.each_col(kvalidx) = solution;

									if (solution.has_nan())
									{
										cout << "NaN here" << "\n";
									}

									break;
								}
							}
						}
						//Place result in main matrix
						if (res.has_nan())
						{
							cout << "why  " << ingam(0) << "  " << sol_k_vec[ingam(0)].size() <<"\n";
						}
						lutcmyk(span(0, 3), span(inds4d(*it), inds4d(*it) + nogrdpts - 1)) = res;
						//res.t().print();
						
						//Mark vertex as solved
						solved(inds3d(*it)) = 1;

						//Fill kmin kmax matrix - all out-of gamut points were solved
						//using one in-gmaut point - refptlstind=ingam(0)
						kminmaxmat(0, inds3d(*it)) = kmins[refptlstind];
						kminmaxmat(1, inds3d(*it)) = kmaxs[refptlstind];
					}
				}
				cout << l << " " << k << " " << j <<"\n";
			} // j
		} // k		
	} //l

	//return lutcmyk;
}

void RevInt::getLabExtmRng(double *labextmarr, double *rangearr)
{
	mat labextm(labextmarr, 3, 2, false, true);
	labextm = m_labextm;
	vec range(rangearr, 3, false, true);
	range = m_labrng;
}

 cmsHTRANSFORM RevInt::makeCForm(const char *profname, const char *luttype)
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

	return hTransform;
}

 cmsHTRANSFORM RevInt::makegammapForm(float *LUTvals, unsigned int nogrdpts)
 {
	 cmsHTRANSFORM hTransformp;
	 cmsHPROFILE hProfilep = cmsCreateProfilePlaceholder(0);

	 cmsSetProfileVersion(hProfilep, 4.2);
	 cmsSetDeviceClass(hProfilep, cmsSigAbstractClass);
	 cmsSetColorSpace(hProfilep, cmsSigRgbData);
	 cmsSetPCS(hProfilep, cmsSigRgbData);

	 cmsFloat32Number* cmsf = LUTvals;
	 cmsStage* cmsStageAllocCLutFloat(cmsContext ContextID,
		 cmsUInt32Number nGridPoints,
		 cmsUInt32Number inputChan,
		 cmsUInt32Number outputChan,
		 const cmsFloat32Number * Table);
	 cmsPipeline* Pipeline = cmsPipelineAlloc(0, 3, 3);
	 cmsPipelineInsertStage(Pipeline, cmsAT_BEGIN, cmsStageAllocCLutFloat(0, nogrdpts, 3, 3, cmsf));
	 cmsWriteTag(hProfilep, cmsSigAToB1Tag, Pipeline);
	 cmsPipelineFree(Pipeline);

	 hTransformp = cmsCreateMultiprofileTransform(&hProfilep,
		 1,
		 TYPE_RGB_FLT,
		 TYPE_RGB_FLT,
		 INTENT_RELATIVE_COLORIMETRIC,
		 0);
	 cmsCloseProfile(hProfilep);

	 return hTransformp;
 }

 cmsHTRANSFORM RevInt::makelab2kForm(float *LUTvals, unsigned int nogrdpts)
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
	 cmsPipeline* Pipeline = cmsPipelineAlloc(0, 3, 2);
	 cmsPipelineInsertStage(Pipeline, cmsAT_BEGIN, cmsStageAllocCLutFloat(0, nogrdpts, 3, 2, cmsf));
	 cmsWriteTag(hProfilep, cmsSigAToB1Tag, Pipeline);
	 cmsPipelineFree(Pipeline);

	 hTransformp = cmsCreateMultiprofileTransform(&hProfilep,
		 1,
		 TYPE_CMYK_DBL,
		 TYPE_CMYK_DBL,
		 INTENT_RELATIVE_COLORIMETRIC,
		 0);
	 cmsCloseProfile(hProfilep);

	 return hTransformp;
 }

 cmsHTRANSFORM RevInt::makelabk2cmykForm(float *LUTvals, size_t nogrdpts)
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

	 return hTransformp;
 }

 void RevInt::applygammapForm(cmsHTRANSFORM gammapform, double* LabIn, double* LabOut, size_t size)
 {
	 cmsDoTransform(gammapform, LabIn, LabOut, size);
 }

 void RevInt::applylab2kForm(cmsHTRANSFORM labk2kform, double* LabIn, double* Kout, size_t size)
 {
	 cmsDoTransform(labk2kform, LabIn, Kout, size);
 }

 void RevInt::applylabk2cmykForm(cmsHTRANSFORM labk2cmykform, double* LabKin, double* CMYKout, size_t size)
 {
	 cmsDoTransform(labk2cmykform, LabKin, CMYKout, size);
 }

 void RevInt::deletegammapForm(cmsHTRANSFORM gammapform)
 {
	 cmsDeleteTransform(gammapform);
 }
 void RevInt::deletelab2kForm(cmsHTRANSFORM labk2kform)
 {
	 cmsDeleteTransform(labk2kform);
 }
 void RevInt::deletelabk2cmykForm(cmsHTRANSFORM labk2cmykform)
 {
	 cmsDeleteTransform(labk2cmykform);
 }

 mat RevInt::permRep(vec x, int nocols, const char* varies)
//Either first or last column varies most rapidly
 {
	 mat A(pow(x.n_elem, nocols), nocols);

	 if (std::string(varies) == "first")
	 {
		 for (int col = 0; col < nocols; col++)
		 {
			 int reps = pow(x.n_elem, col);
			 int row = 0;
			 while (row < A.n_rows)
			 {
				 for (int el = 0; el < x.n_elem; el++)
				 {
					 for (int rep = 0; rep < reps; rep++)
					 {
						 A(row, col) = x(el);
						 row++;
					 }
				 }
			 }
		 }
	 }
	 else if (std::string(varies) == "last")
	 {
		 for (int col = nocols-1; col >= 0; col--)
		 {
			 int reps = pow(x.n_elem, col);
			 int row = 0;
			 while (row < A.n_rows)
			 {
				 for (int el = 0; el < x.n_elem; el++)
				 {
					 for (int rep = 0; rep < reps; rep++)
					 {
						 A(row, col) = x(el);
						 row++;
					 }
				 }
			 }
		 }
	 }

	 return A;
 }

 Mat<int> RevInt::permWoRep(int x[], int size)
 {
	 Mat<int>A(factorial(size), size);
	 int row = 0;
	 do {
		 A.row(row) = s32_rowvec(x, size);
		 row++;
	 } while (std::next_permutation(x, x + 4));
	 return A;
 }

 int RevInt::factorial(int n)
 {
	 return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
 }

 int RevInt::lutInd3D(int i, int j, int k, int nopts)
 {
	 int ind = pow(nopts,2)*k + nopts*j + i;
	 return ind;
 }

 int RevInt::lutInd4D(int i, int j, int k, int l, int nopts)
 {
	 int ind = pow(nopts, 3)*l + pow(nopts, 2)*k + nopts * j + i;
	 return ind;
 }

 Mat<int> RevInt::intersect(Mat<int> A, Mat<int> B)
 {
	 Mat<int> C;

	 for (int rowA = 0; rowA < A.n_rows; rowA++)
	 {
		 bool match = false;
		 int rowB = 0;
		 while ((rowB < B.n_rows) && match == false)
		 {
			 if (arma::approx_equal(A.row(rowA), B.row(rowB), "absdiff", 0.1))
			 {
				 C.resize(C.n_rows+1, A.n_cols);
				 C.row(C.n_rows-1) = A.row(rowA);
				 match = true;
			 }				 
			 rowB++;
		 }
	 }
	 return C;
 }

 vec RevInt::overlap(vec c0, vec c1)
 {
	 double cmin = 0;
	 double cmax = 0;
	 vec interval;

	 double min1 = std::min(c0(0), c1(0));
	 double max1 = std::max(c0(0), c1(0));
	 double min2 = std::min(c0(1), c1(1));
	 double max2 = std::max(c0(1), c1(1));

	 if (min2 > max1) 
	 {
		 return interval;
	 }
	 else if (max2 < min1)
	 {
		 return interval;
	 }
	 else
	 {
		 cmin = std::max(min1, min2);
		 cmax = std::min(max1, max2);
		 interval.set_size(2);
		 interval(0) = cmin;
		 interval(1) = cmax;
	 }

	 if (c0.n_elem>=3)
	 {
		 for (int i = 2; i < c0.n_elem; i++)
		 {
			 min1 = cmin;
			 max1 = cmax;
			 cmin = 0;
			 cmax = 0;
			 interval.reset();
			 min2 = std::min(c0(i), c1(i));
			 max2 = std::max(c0(i), c1(i));

			 if (min2 > max1)
			 {
				 return interval;
			 }
			 else if (max2 < min1)
			 {
				 return interval;
			 }
			 else
			 {
				 interval.set_size(2);
				 cmin = std::max(min1, min2);
				 cmax = std::min(max1, max2);
				 interval(0) = cmin;
				 interval(1) = cmax;
			 }		 
		 }
	 }
	 return interval;
 }

RevInt::~RevInt()
{
	cmsDeleteTransform(m_cfAtoB);
	cmsDeleteTransform(m_cfBtoA);
}
