#include <mex.h> 
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <utility>

#include <dirent.h>

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/dirname.h>

#include <igl/readMESH.h>
#include <igl/readOBJ.h>

#include <Eigen/Core>

using namespace Eigen;
using namespace std;

void read_mesh_files_from_directory(const std::string& directory, std::vector<std::string>& mesh_list){

	using namespace std;

	std::vector<std::string> file_list;

	DIR *d;
	struct dirent *dir;
	int i=0;
	d = opendir(directory.c_str());
	if (d)
	{
		while ((dir = readdir(d)) != NULL)
		{
			i++;
			file_list.push_back(dir->d_name);
		}
		closedir(d);
	}
	for(int i=0; i<file_list.size(); i++){
		std::string file = file_list[i];
		size_t last_index = file.find_last_of(".");
		if(file.substr(last_index+1) == "mesh"){
			std::string name = file.substr(0,last_index);
			mesh_list.push_back(name);
		}
	}
}



void read_blendshape_list(const std::string& directory, const std::string& list_file,
		std::vector<std::string>& name_list,
		std::vector<Eigen::MatrixXd>& V_list)
{
	FILE* file;
	file= fopen(list_file.c_str(), "rb");
	int num_list;
	fscanf(file, "%d\n", &num_list);


	std::vector<int> reordered_id(num_list);
	for(int i=0; i<num_list; i++)
	{
		//std::string bs_name;
		char bs_char[256];
		fscanf(file, "%s\n", bs_char);
		std::string bs_name(bs_char);
		for(int j=0; j<name_list.size(); j++){
			if(bs_name==name_list[j]){
				reordered_id[i]=j;
			}
		}
	}
	std::vector<std::string> reorder_name_list(num_list);
	for(int i=0; i<reordered_id.size(); i++){
		reorder_name_list[i] = name_list[reordered_id[i]];
	}
	name_list = reorder_name_list;
	fclose(file);

	V_list.resize(name_list.size());
	for(int j=0; j<name_list.size(); j++){
		// std::cout<<"name_list["<<j<<"] : "<<name_list[j]<<std::endl;
		std::string mesh_path = directory+"/"+name_list[j]+".mesh";
		Eigen::MatrixXd BV;
		Eigen::MatrixXi BT, BF;
		igl::readMESH(mesh_path, BV,BT,BF);
		// igl::readOBJ(mesh_path, BV,BF);
		V_list[j] = BV;
	}
}




void read_blendshape_anim(const std::string& anim_file, std::vector<Eigen::VectorXd>& w_list)
{
	FILE* file;
	file= fopen(anim_file.c_str(), "rb");
	int num_bs, num_frame;
	fscanf(file, "%d %d\n", &num_bs, &num_frame);

	w_list.resize(num_frame);
	for(int i=0; i<num_frame; i++){
		Eigen::VectorXd w(num_bs);
		for(int j=0; j<num_bs; j++){
			double val;
			fscanf(file, "%lf", &val);
			w(j) = val;
		}
		w_list[i] = w;
	}
	fclose(file);
}



void construct_blendshape_mat(const Eigen::MatrixXd& V,
		const std::vector<Eigen::MatrixXd>& BV_list,
		//Eigen::RowVectorXd Vcenter, double Vbound,
		Eigen::SparseMatrix<double>& BS)
{
	Eigen::MatrixXd VT = V.transpose();
	Eigen::VectorXd VCol(Eigen::Map<Eigen::VectorXd>(VT.data(), V.cols()*V.rows()));
	std::vector<Eigen::Triplet<double> > coefficients;
	for(int j=0; j<BV_list.size(); j++){
		Eigen::MatrixXd BV = BV_list[j];

		//normalize_obj(BV, Vcenter, Vbound);

		Eigen::MatrixXd BVT = BV.transpose();

		Eigen::VectorXd bVcol = Eigen::Map<Eigen::VectorXd>(BVT.data(), BV.cols()*BV.rows());

		Eigen::VectorXd deltaV = bVcol-VCol;
		for(int i=0; i<deltaV.size(); i++){
			if(deltaV(i)>1e-12){
				coefficients.push_back(Eigen::Triplet<double>(i,j, deltaV(i)));
			}
		}
	}
	BS.resize(VCol.size(), BV_list.size());
	BS.setFromTriplets(coefficients.begin(), coefficients.end());
}

void construct_blendshape_mat(const Eigen::MatrixXd& V,
		const std::vector<Eigen::MatrixXd>& BV_list,
		//Eigen::RowVectorXd Vcenter, double Vbound,
		Eigen::MatrixXd& BS)
{
	Eigen::MatrixXd VT = V.transpose();
	Eigen::VectorXd VCol(Eigen::Map<Eigen::VectorXd>(VT.data(), V.cols()*V.rows()));

	BS.resize(VCol.size(), BV_list.size());
	BS.setZero();
	for(int j=0; j<BV_list.size(); j++){
		Eigen::MatrixXd BV = BV_list[j];

		//normalize_obj(BV, Vcenter, Vbound);

		Eigen::MatrixXd BVT = BV.transpose();

		Eigen::VectorXd bVcol = Eigen::Map<Eigen::VectorXd>(BVT.data(), BV.cols()*BV.rows());

		Eigen::VectorXd deltaV = bVcol-VCol;
		for(int i=0; i<deltaV.size(); i++){
			//if(deltaV(i)>1e-12){
				BS(i,j)=deltaV(i);
			//}
		}
	}
}

void read_blendshape_data(
	std::string bs_list_file, std::string bs_folder, std::string animfile,
	const Eigen::MatrixXd& TV,
	//Eigen::SparseMatrix<double>& BS,
	Eigen::MatrixXd& BS,
	std::vector<Eigen::VectorXd>& T_list
	)
{

	std::vector<std::string> mesh_name_list;
	read_mesh_files_from_directory(bs_folder, mesh_name_list);

	std::vector<Eigen::MatrixXd> BV_list;
	read_blendshape_list(bs_folder, bs_list_file, mesh_name_list, BV_list);
	read_blendshape_anim(animfile, T_list);

	//Eigen::SparseMatrix<double> BS;
	//construct_blendshape_mat(TV, BV_list, Vcenter, Vbound, BS);
	construct_blendshape_mat(TV, BV_list, BS);
}

void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[])
{
    Eigen::MatrixXd TV;

    char* bs_list_file = mxArrayToString(prhs[0]);
    char* bs_folder = mxArrayToString(prhs[1]);
    char* animfile = mxArrayToString(prhs[2]);
	igl::matlab::parse_rhs_double(prhs+3,TV);
	//Eigen::SparseMatrix<double> BS;
	Eigen::MatrixXd BS;
    std::vector<Eigen::VectorXd> T_list;


	read_blendshape_data(bs_list_file, bs_folder, animfile, TV, BS, T_list);

	igl::matlab::prepare_lhs_double(BS, plhs);

    plhs[1] = mxCreateCellMatrix(T_list.size(),1);

    mxArray *x;
    for(int i=0; i<T_list.size(); i++){
        const int m = T_list[i].rows();
        const int n = T_list[i].cols();
        x = mxCreateDoubleMatrix(m,n, mxREAL);
        Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > map(mxGetPr(x),m,n);
        map = T_list[i].template cast<double>();
        mxSetCell(plhs[1], i, x);
    }
    return;
}



