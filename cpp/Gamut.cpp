#include "Gamut.h"


Gamut::Gamut(double *Lab, size_t no_Lab, int *ch, size_t no_ch)
{
	m_no_Lab = no_Lab;
	m_no_ch = no_ch;

	mat_Lab = mat(Lab, 3, m_no_Lab);
	mat_ch = Mat<int>(ch, 3, m_no_ch);
	uword Imax = mat_Lab.row(0).index_max();
	uword Imin = mat_Lab.row(0).index_min();
	wpt = mat_Lab.submat(0, Imax, 2, Imax);
	bpt = mat_Lab.submat(0, Imin, 2, Imin);
	centroid = (wpt + bpt) / 2;

	uvec vec_ch = conv_to<uvec>::from(arma::vectorise(mat_ch));
	mat positions = mat_Lab.cols(vec_ch);
	positions.each_col() -= centroid;

	//Map of vectors setup
	vec_map.insert(std::make_pair(1, oct1));
	vec_map.insert(std::make_pair(2, oct2));
	vec_map.insert(std::make_pair(3, oct3));
	vec_map.insert(std::make_pair(4, oct4));
	vec_map.insert(std::make_pair(5, oct5));
	vec_map.insert(std::make_pair(6, oct6));
	vec_map.insert(std::make_pair(7, oct7));
	vec_map.insert(std::make_pair(8, oct8));

	int face_no = 0;
	for (int j = 0; j < positions.n_cols; j+=3)
	{
		//Get positions (octants) for three points of single hull face
		int first_position = 0;
		int second_position = 0;
		int third_position = 0;
		first_position = GetOctand(positions.col(j));
		second_position = GetOctand(positions.col(j+1));
		third_position = GetOctand(positions.col(j+2));
		//Add face to corresponding vector based on first point position
		vec_map[first_position].push_back(face_no);
		//Add face to corresponding vector based on second and third point positions
		//Skip if their position equals previously added position.
		if (second_position != first_position)
			vec_map[second_position].push_back(face_no);
		if (third_position != first_position && third_position != second_position)
			vec_map[third_position].push_back(face_no);
		face_no++;
	}

	std::ofstream out("GCR_py_wrap_Gamut_constructor.txt");

	out << "mat_Lab: " << mat_Lab.n_rows << " x " << mat_Lab.n_cols << "\n" << 
		"no_Lab: " << no_Lab << "\n" <<
		"mat_ch: " << mat_ch.n_rows << " x " << mat_ch.n_cols << "\n" <<
		"no_ch: " << no_ch << "\n" <<
		"centroid: " << centroid.n_rows << " x " << centroid.n_cols << "\n";
	out.close();
}

void Gamut::map(double *Lab_set, double *Lab_set_new, size_t no_Lab_set, unsigned short *outgamarr, pfnPrintProgress PrintCallback)
{
	mat mat_Lab_set(Lab_set, 3, no_Lab_set, false);
	mat mat_Lab_set_new(Lab_set_new, 3, no_Lab_set, false);
	mat_Lab_set_new = mat_Lab_set;

	//Set outgamarr elements to zero
	for (size_t i = 0; i < no_Lab_set; i++)
	{
		*(outgamarr + i) = 0;
	}

	//Get positions and octants of colors
	mat mat_Lab_set_positions = mat_Lab_set;
	mat_Lab_set_positions.each_col() -= centroid;
	std::vector<int> octants(mat_Lab_set.n_cols);
	for (size_t i = 0; i < mat_Lab_set.n_cols; i++)
	{
		octants[i] = GetOctand(mat_Lab_set_positions.col(i));
	}

	//For each color in mat_Lab_set
	for (uword i = 0; i < mat_Lab_set.n_cols; i++)
	{
		vec pta;
		vec ptac;
		try
		{
			pta = mat_Lab_set.col(i);
			ptac = pta - centroid;
		}
		catch (const std::exception& ex)
		{
			std::ofstream out("GCR_py_wrap_Gamut33.txt");
			std::string str(ex.what());
			out << str << "\n" << 
				"mat_Lab_set: " << mat_Lab_set.n_rows << " x " << mat_Lab_set.n_cols << "\n" << 
				"no_Lab_set: " << no_Lab_set << "\n" <<
				"centroid: " << centroid.n_rows << " x " << centroid.n_cols << "\n" <<
				"mat_Lab: " << mat_Lab.n_rows << " x " << mat_Lab.n_cols << "\n" <<
				"mat_ch: " << mat_ch.n_rows << " x " << mat_ch.n_cols << "\n" <<
				"m_no_Lab: " << m_no_Lab << "\n" <<
				"m_no_ch: " << m_no_ch << "\n";
			out.close();
		}
		

		mat sol(3, mat_ch.n_cols);
		uvec intersections(mat_ch.n_cols, fill::zeros);
		uvec outgam(mat_ch.n_cols, fill::zeros);

		//Seek solution in current color's octant
		for (uword j = 0; j < vec_map[octants[i]].size(); j++)
		{
			mat A(3, 3);
			
			try
			{
				A.col(0) = -ptac;
				A.col(1) = mat_Lab.col(mat_ch(1, vec_map[octants[i]][j])) - mat_Lab.col(mat_ch(0, vec_map[octants[i]][j]));
				A.col(2) = mat_Lab.col(mat_ch(2, vec_map[octants[i]][j])) - mat_Lab.col(mat_ch(0, vec_map[octants[i]][j]));
				sol.col(j) = arma::solve(A, centroid - mat_Lab.col(mat_ch(0, vec_map[octants[i]][j])));
			}
			catch (const std::exception& ex)
			{
				std::ofstream out("GCR_py_wrap_Gamut55.txt");
				std::string str(ex.what());
				out << str;
				out.close();
			}
			

			if (sol(1, j) >= 0 && sol(1, j) <= 1 && sol(2, j) >= 0 && sol(2, j) <= 1 && sol(1, j) + sol(2, j) <=1)
			{
				intersections(j) = 1;
			}

			if (sol(0, j) > 0 && sol(0, j) < 1)
			{
				outgam(j) = 1;				
			}						
		}
		uvec intersetction = find((intersections % outgam) == 1);

		if (intersetction.n_elem > 1)
		{
			intersetction.resize(1);
		}

		if (intersetction.is_empty() == false)
		{
			mat_Lab_set_new.col(i) = (ptac * sol.submat(zeros<uvec>(1), intersetction).eval()(0, 0)) + centroid;
			*(outgamarr + i) = 1;
			//cout << "\n";
			//mat_Lab_set.col(i).print();
			//mat_Lab_set_new.col(i).print();
		}
		/*if (std::fmod((double(i)/mat_Lab_set.n_cols*100.0), 5) == 0)
		{
			PrintCallback(size_t(i), size_t(mat_Lab_set.n_cols));
		}*/
		PrintCallback(size_t(i), size_t(mat_Lab_set.n_cols));
	}
}

void Gamut::PrintProgress(pfnPrintProgress PrintCallback, size_t current, size_t total)
{
	PrintCallback(current, total);
}

int Gamut::GetOctand(vec Lab_pos)
{
	if (Lab_pos(0) < 0 && Lab_pos(1) < 0 && Lab_pos(2) < 0)
		return 1;
	else if (Lab_pos(0) < 0 && Lab_pos(1) < 0 && Lab_pos(2) > 0)
		return 2;
	else if (Lab_pos(0) < 0 && Lab_pos(1) > 0 && Lab_pos(2) < 0)
		return 3;
	else if (Lab_pos(0) < 0 && Lab_pos(1) > 0 && Lab_pos(2) > 0)
		return 4;
	else if (Lab_pos(0) > 0 && Lab_pos(1) < 0 && Lab_pos(2) < 0)
		return 5;
	else if (Lab_pos(0) > 0 && Lab_pos(1) < 0 && Lab_pos(2) > 0)
		return 6;
	else if (Lab_pos(0) > 0 && Lab_pos(1) > 0 && Lab_pos(2) < 0)
		return 7;
	else
		return 8;
}

Gamut::~Gamut()
{

}
