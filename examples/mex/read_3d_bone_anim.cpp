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
}


void read_bone_anim(std::string anim_file,
    const Eigen::MatrixXd& C,
    const Eigen::MatrixXi& BE,
    const Eigen::VectorXi& P,
    std::vector<Eigen::MatrixXd>& T_list)
{
  FILE* file;
  file= fopen(anim_file.c_str(), "rb");

  double degree2radian = igl::PI/180.0;

  int num_bone, num_frame;
  double val = 0;

  typedef std::vector<Eigen::Quaterniond,
      Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


  fscanf(file, "%d %d\n", &num_bone, &num_frame);
  T_list.resize(num_frame);
  RotationList rot_list(num_bone);
  std::vector<Eigen::Vector3d> tran_list(num_bone);
  std::vector<Eigen::Affine3d> rest_list(num_bone);

  Eigen::Vector3d root_rest_tran;
  Eigen::Affine3d root_affine;
  for (int j = 0; j < 3; j++) {
    fscanf(file, "%lf", &val);
    root_rest_tran(j) = val;
  }
  Eigen::Vector3d euler;
  for (int j = 0; j < 3; j++) {
    fscanf(file, "%lf", &val);
    euler(j) = val;
  }
  euler = euler.array() * degree2radian;
  Eigen::Quaterniond q;
  euler_to_quat(euler, q);
  root_affine = Eigen::Affine3d::Identity();
  root_affine.rotate(q);

  for (int i = 0; i < num_bone; i++) {
    Eigen::Vector3d euler;
    for (int j = 0; j < 3; j++) {
      fscanf(file, "%lf", &val);
      euler(j) = val;
    }
    euler = euler.array() * degree2radian;
    Eigen::Quaterniond q;
    euler_to_quat(euler, q);
    Eigen::Affine3d a = Eigen::Affine3d::Identity();
    a.rotate(q);
    rest_list[i] = a;
  }
  //Eigen::Affine3d root_affine = rest_list[0];
  Eigen::Vector3d root_tran, root_rot;
  for (int k = 0; k < num_frame; k++) {
    for (int j = 0; j < 3; j++) {
      fscanf(file, "%lf", &val);
      root_rot(j) = val;
    }
    for (int j = 0; j < 3; j++) {
      fscanf(file, "%lf", &val);
      root_tran(j) = val;
    }
    root_rot = root_rot.array() * degree2radian;
    Eigen::Quaterniond root_q;
    euler_to_quat(root_rot, root_affine, root_q);

    for (int i = 0; i < num_bone; i++) {
      Eigen::Vector3d euler;
      for (int j = 0; j < 3; j++) {
        fscanf(file, "%lf", &val);
        euler(j) = val;
      }
      euler = euler.array() * degree2radian;
      Eigen::Quaterniond q;
      euler_to_quat(euler, rest_list[i], q);

      if(P(i)==-1){
        int root_cnt = 0;
        for(int ii=0; ii<num_bone; ii++)
          if(P(ii)==-1)
            root_cnt++;
        if(root_cnt>1)
          q = q*root_q;
        tran_list[i] = root_tran - root_rest_tran;
      }
      else
        tran_list[i] = Eigen::Vector3d::Zero();
      rot_list[i] = q;
    }
    RotationList vQ;
    std::vector<Eigen::Vector3d> vT;
    igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

    Eigen::MatrixXd T(num_bone * 4, 3);
    for (int i = 0; i < num_bone; i++) {
      Eigen::Affine3d a = Eigen::Affine3d::Identity();
      a.translate(vT[i]);
      a.rotate(vQ[i]);
      T.block(i * 4, 0, 4, 3) = a.matrix().transpose().block(0, 0, 4, 3);
    }
    T_list[k] = T;
  }
  fclose(file);
}

void read_bone_anim(std::string anim_file,
    const Eigen::MatrixXd& C,
    const Eigen::MatrixXi& BE,
    const Eigen::VectorXi& P,
    const Eigen::VectorXd& center,
    const double& scale,
    std::vector<Eigen::MatrixXd>& T_list)
{
  FILE* file;
  file= fopen(anim_file.c_str(), "rb");

  double degree2radian = igl::PI/180.0;

  int num_bone, num_frame;
  double val = 0;

  typedef std::vector<Eigen::Quaterniond,
      Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


  fscanf(file, "%d %d\n", &num_bone, &num_frame);
  T_list.resize(num_frame);
  RotationList rot_list(num_bone);
  std::vector<Eigen::Vector3d> tran_list(num_bone);
  std::vector<Eigen::Affine3d> rest_list(num_bone);

  Eigen::Vector3d root_rest_tran;
  Eigen::Affine3d root_affine;
  for (int j = 0; j < 3; j++) {
    fscanf(file, "%lf", &val);
    root_rest_tran(j) = val;
  }
  root_rest_tran = (root_rest_tran-center)/scale;

  Eigen::Vector3d euler;
  for (int j = 0; j < 3; j++) {
    fscanf(file, "%lf", &val);
    euler(j) = val;
  }
  euler = euler.array() * degree2radian;
  Eigen::Quaterniond q;
  euler_to_quat(euler, q);
  root_affine = Eigen::Affine3d::Identity();
  root_affine.rotate(q);

  for (int i = 0; i < num_bone; i++) {
    Eigen::Vector3d euler;
    for (int j = 0; j < 3; j++) {
      fscanf(file, "%lf", &val);
      euler(j) = val;
    }
    euler = euler.array() * degree2radian;
    Eigen::Quaterniond q;
    euler_to_quat(euler, q);
    Eigen::Affine3d a = Eigen::Affine3d::Identity();
    a.rotate(q);
    rest_list[i] = a;
  }
  //Eigen::Affine3d root_affine = rest_list[0];
  Eigen::Vector3d root_tran, root_rot;
  for (int k = 0; k < num_frame; k++) {
    for (int j = 0; j < 3; j++) {
      fscanf(file, "%lf", &val);
      root_rot(j) = val;
    }
    for (int j = 0; j < 3; j++) {
      fscanf(file, "%lf", &val);
      root_tran(j) = val;
    }
    root_tran = (root_tran - center)/scale;
    root_rot = root_rot.array() * degree2radian;
    Eigen::Quaterniond root_q;
    euler_to_quat(root_rot, root_affine, root_q);

    for (int i = 0; i < num_bone; i++) {
      Eigen::Vector3d euler;
      for (int j = 0; j < 3; j++) {
        fscanf(file, "%lf", &val);
        euler(j) = val;
      }
      euler = euler.array() * degree2radian;
      Eigen::Quaterniond q;
      euler_to_quat(euler, rest_list[i], q);

      if(P(i)==-1){
        int root_cnt = 0;
        for(int ii=0; ii<num_bone; ii++)
          if(P(ii)==-1)
            root_cnt++;
        if(root_cnt>1)
          q = q*root_q;
        tran_list[i] = root_tran - root_rest_tran;
      }
      else
        tran_list[i] = Eigen::Vector3d::Zero();
      rot_list[i] = q;
    }
    RotationList vQ;
    std::vector<Eigen::Vector3d> vT;
    igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

    Eigen::MatrixXd T(num_bone * 4, 3);
    for (int i = 0; i < num_bone; i++) {
      Eigen::Affine3d a = Eigen::Affine3d::Identity();
      a.translate(vT[i]);
      a.rotate(vQ[i]);
      T.block(i * 4, 0, 4, 3) = a.matrix().transpose().block(0, 0, 4, 3);
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
    igl::matlab::parse_rhs_index(prhs+2, BE);

    igl::matlab::parse_rhs_double(prhs+3, center);
    scale = (double) *mxGetPr(prhs[4]);

    igl::directed_edge_parents(BE, P);

    std::vector<Eigen::MatrixXd> T_list;
    read_bone_anim(anim_file, C,BE,P,center,scale, T_list);
    //read_bone_anim(anim_file, C,BE,P, T_list);


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