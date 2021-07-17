#include <mex.h> 
#include <math.h>
#include <iostream>
#include <vector>

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>

#include <igl/PI.h>
#include <igl/forward_kinematics.h>
#include <igl/directed_edge_parents.h>

#include <Eigen/Core>

using namespace Eigen;
using namespace std;

void euler_to_quat(const Eigen::Vector3d& euler, Eigen::Quaterniond& q) {

  q = Eigen::AngleAxisd(euler(2), Eigen::Vector3d::UnitZ())
      * Eigen::AngleAxisd(euler(1), Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd(euler(0), Eigen::Vector3d::UnitX());
}


void euler_to_quat(const Eigen::Vector3d& euler, const Eigen::Affine3d& a, Eigen::Quaterniond& q) {
  q = Eigen::AngleAxisd(euler(2), a.rotation().col(2))
      * Eigen::AngleAxisd(euler(1), a.rotation().col(1))
      * Eigen::AngleAxisd(euler(0), a.rotation().col(0));
  // q = Eigen::AngleAxisd(euler(0), a.rotation().col(0))
  //     * Eigen::AngleAxisd(euler(1), a.rotation().col(1))
  //     * Eigen::AngleAxisd(euler(2), a.rotation().col(2));
}


void read_pnt_anim(std::string anim_file,
    const Eigen::MatrixXd& C,
    const Eigen::VectorXd& center,
    const double& scale,
    std::vector<Eigen::MatrixXd>& T_list)
{
  FILE* file;
  file= fopen(anim_file.c_str(), "rb");

  int dim = 3;//= C.cols();

  double degree2radian = igl::PI/180.0;

  int num_bone, num_frame;
  double val = 0;

  typedef std::vector<Eigen::Quaterniond,
      Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


  fscanf(file, "%d %d\n", &num_bone, &num_frame);
  T_list.resize(num_frame);
  //std::cout<<"num_bone: "<<num_bone<<", num_frame: "<<num_frame<<std::endl;

  std::vector<Eigen::Vector3d> rest_tran_list(num_bone);
  std::vector<Eigen::Affine3d> rest_affine_list(num_bone);

  for(int i=0; i<num_bone ; i++){
    Eigen::Vector3d rest_tran;
    for (int j = 0; j < dim; j++) {
      fscanf(file, "%lf", &val);
      rest_tran(j) = val;
    }
    rest_tran_list[i] = (rest_tran-center)/scale;

    Eigen::Vector3d euler;
    for (int j = 0; j < dim; j++) {
      fscanf(file, "%lf", &val);
      euler(j) = val;
    }
    euler = euler.array() * degree2radian;
    Eigen::Quaterniond q;
    euler_to_quat(euler, q);
    Eigen::Affine3d root_affine = Eigen::Affine3d::Identity();
    root_affine.rotate(q);
    rest_affine_list[i] = root_affine;
    //std::cout<<"i: "<<i<<", mat: "<<root_affine.matrix()<<std::endl;
  }
  RotationList rot_list(num_bone);
  std::vector<Eigen::Vector3d> tran_list(num_bone);

  for (int k = 0; k < num_frame; k++) {
    for (int i = 0; i < num_bone; i++) {
      Eigen::Vector3d tran;
      for (int j = 0; j < dim; j++) {
        fscanf(file, "%lf", &val);
        tran(j) = val;
      }
      tran = (tran - center)/scale;
      Eigen::Vector3d euler;
      for (int j = 0; j < dim; j++) {
        fscanf(file, "%lf", &val);
        euler(j) = val;
      }
      euler = euler.array() * degree2radian;
      Eigen::Quaterniond q;
      euler_to_quat(euler, rest_affine_list[i], q);
      //std::cout<<"1 i: "<<i<<", mat: "<<rest_affine_list[i].rotation()<<std::endl;

      q = rest_affine_list[i].rotation().transpose()*q;

      //euler_to_quat(euler, q);

      //std::cout<<"q w: "<<q.w()<<"< vec: "<<q.vec()<<std::endl;

      const Eigen::Vector3d cn = C.row(i).transpose();

      tran_list[i] = tran - rest_tran_list[i] + cn-q*cn;
      //tran_list[i] = tran - rest_tran_list[i];
      rot_list[i] = q;
    }
    RotationList vQ = rot_list;
    std::vector<Eigen::Vector3d> vT= tran_list;
    //igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

    Eigen::MatrixXd T(num_bone * (dim+1), dim);
    for (int i = 0; i < num_bone; i++) {
      Eigen::Affine3d a = Eigen::Affine3d::Identity();
      a.translate(vT[i]);
      a.rotate(vQ[i]);
      T.block(i * (dim+1), 0, dim+1, dim) = a.matrix().transpose().block(0, 0, dim+1, dim);
    }
    T_list[k] = T;
  }
  fclose(file);

}


// void read_pnt_anim(std::string anim_file,
//     std::vector<Eigen::MatrixXd>& T_list)
// {
//   FILE* file;
//   file= fopen(anim_file.c_str(), "rb");

//   int dim = 3;//= C.cols();

//   double degree2radian = igl::PI/180.0;

//   int num_bone, num_frame;
//   double val = 0;

//   typedef std::vector<Eigen::Quaterniond,
//       Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


//   fscanf(file, "%d %d\n", &num_bone, &num_frame);
//   T_list.resize(num_frame);
//   //std::cout<<"num_bone: "<<num_bone<<", num_frame: "<<num_frame<<std::endl;

//   std::vector<Eigen::Vector3d> rest_tran_list(num_bone);
//   std::vector<Eigen::Affine3d> rest_affine_list(num_bone);

//   for(int i=0; i<num_bone ; i++){
//     Eigen::Vector3d rest_tran;
//     for (int j = 0; j < dim; j++) {
//       fscanf(file, "%lf", &val);
//       rest_tran(j) = val;
//     }
//     rest_tran_list[i] = rest_tran;

//     Eigen::Vector3d euler;
//     for (int j = 0; j < dim; j++) {
//       fscanf(file, "%lf", &val);
//       euler(j) = val;
//     }
//     euler = euler.array() * degree2radian;
//     Eigen::Quaterniond q;
//     euler_to_quat(euler, q);
//     Eigen::Affine3d root_affine = Eigen::Affine3d::Identity();
//     root_affine.rotate(q);
//     rest_affine_list[i] = root_affine;
//   }
//   RotationList rot_list(num_bone);
//   std::vector<Eigen::Vector3d> tran_list(num_bone);

//   for (int k = 0; k < num_frame; k++) {
//     for (int i = 0; i < num_bone; i++) {
//       Eigen::Vector3d tran;
//       for (int j = 0; j < dim; j++) {
//         fscanf(file, "%lf", &val);
//         tran(j) = val;
//       }
//       Eigen::Vector3d euler;
//       for (int j = 0; j < dim; j++) {
//         fscanf(file, "%lf", &val);
//         euler(j) = val;
//       }
//       euler = euler.array() * degree2radian;
//       Eigen::Quaterniond q;
//       euler_to_quat(euler, rest_affine_list[i], q);

//       tran_list[i] = tran - rest_tran_list[i];
//       rot_list[i] = q;
//     }
//     RotationList vQ = rot_list;
//     std::vector<Eigen::Vector3d> vT= tran_list;
//     //igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

//     Eigen::MatrixXd T(num_bone * (dim+1), dim);
//     for (int i = 0; i < num_bone; i++) {
//       Eigen::Affine3d a = Eigen::Affine3d::Identity();
//       a.translate(vT[i]);
//       a.rotate(vQ[i]);
//       T.block(i * (dim+1), 0, dim+1, dim) = a.matrix().transpose().block(0, 0, dim+1, dim);
//     }
//     T_list[k] = T;
//   }
//   fclose(file);
// }

void read_pnt_anim(std::string anim_file,
    Eigen::MatrixXd C,
    std::vector<Eigen::MatrixXd>& T_list)
{
  FILE* file;
  file= fopen(anim_file.c_str(), "rb");

  int dim = 3;//= C.cols();

  double degree2radian = igl::PI/180.0;

  int num_bone, num_frame;
  double val = 0;

  typedef std::vector<Eigen::Quaterniond,
      Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


  fscanf(file, "%d %d\n", &num_bone, &num_frame);
  T_list.resize(num_frame);
  //std::cout<<"num_bone: "<<num_bone<<", num_frame: "<<num_frame<<std::endl;

  std::vector<Eigen::Vector3d> rest_tran_list(num_bone);
  std::vector<Eigen::Affine3d> rest_affine_list(num_bone);

  for(int i=0; i<num_bone ; i++){
    Eigen::Vector3d rest_tran;
    for (int j = 0; j < dim; j++) {
      fscanf(file, "%lf", &val);
      rest_tran(j) = val;
    }
    rest_tran_list[i] = rest_tran;

    Eigen::Vector3d euler;
    for (int j = 0; j < dim; j++) {
      fscanf(file, "%lf", &val);
      euler(j) = val;
    }
    euler = euler.array() * degree2radian;
    Eigen::Quaterniond q;
    euler_to_quat(euler, q);
    Eigen::Affine3d root_affine = Eigen::Affine3d::Identity();
    root_affine.rotate(q);
    rest_affine_list[i] = root_affine;
  }
  RotationList rot_list(num_bone);
  std::vector<Eigen::Vector3d> tran_list(num_bone);

  for (int k = 0; k < num_frame; k++) {
    for (int i = 0; i < num_bone; i++) {
      Eigen::Vector3d tran;
      for (int j = 0; j < dim; j++) {
        fscanf(file, "%lf", &val);
        tran(j) = val;
      }
      Eigen::Vector3d euler;
      for (int j = 0; j < dim; j++) {
        fscanf(file, "%lf", &val);
        euler(j) = val;
      }
      euler = euler.array() * degree2radian;
      Eigen::Quaterniond q;
      euler_to_quat(euler, rest_affine_list[i], q);

      q = rest_affine_list[i].rotation().transpose()*q;

      const Eigen::Vector3d cn = C.row(i).transpose();

      tran_list[i] = tran - rest_tran_list[i] + cn-q*cn;
      rot_list[i] = q;//cn-q*cn;//q;
    }
    RotationList vQ = rot_list;
    std::vector<Eigen::Vector3d> vT= tran_list;
    //igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

    Eigen::MatrixXd T(num_bone * (dim+1), dim);
    for (int i = 0; i < num_bone; i++) {
      Eigen::Affine3d a = Eigen::Affine3d::Identity();
      a.translate(vT[i]);
      a.rotate(vQ[i]);
      T.block(i * (dim+1), 0, dim+1, dim) = a.matrix().transpose().block(0, 0, dim+1, dim);
    }
    T_list[k] = T;
  }
  fclose(file);
}


void mexFunction(
  int          nlhs,
  mxArray      *plhs[],
  int          nrhs,
  const mxArray *prhs[])
{
    Eigen::MatrixXd C;
    Eigen::MatrixXi BE;
    Eigen::VectorXi P;
    Eigen::VectorXd center;
    double scale;
    char* anim_file = mxArrayToString(prhs[0]);

    igl::matlab::parse_rhs_double(prhs+1, C);
    igl::matlab::parse_rhs_double(prhs+2, center);
    scale = (double) *mxGetPr(prhs[3]);




    std::vector<Eigen::MatrixXd> T_list;
    //read_pnt_anim(anim_file, center, scale, T_list);
    //read_pnt_anim(anim_file, C, T_list);
    read_pnt_anim(anim_file, C, center, scale, T_list);


    plhs[0] = mxCreateCellMatrix(T_list.size(),1);

    mxArray *x;
    for(int i=0; i<T_list.size(); i++){
        const int m = T_list[i].rows();
        const int n = T_list[i].cols();
        x = mxCreateDoubleMatrix(m,n, mxREAL);
        Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > map(mxGetPr(x),m,n);
        map = T_list[i].template cast<double>();
        mxSetCell(plhs[0], i, x);
    }
    return;
}