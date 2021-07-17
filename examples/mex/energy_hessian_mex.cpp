#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;


double amips_param;
double c1_param;
double c2_param;
double d1_param;


Vector2d energy_derivative(int type, Vector2d S)
{
    
    Vector2d es;
    double J;
    double l1;
    double l2;
    
    switch(type)
    {
        case 0: //arap
            es[0] = 2 * (S[0] - 1);
            es[1] = 2 * (S[1] - 1);
            break;
            
        case 1: // mips
            es[0] = 1.0 / S[1] - S[1] / (S[0] * S[0]);
            es[1] = 1.0 / S[0] - S[0] / (S[1] * S[1]);
            break;
            
        case 2: // iso
            es[0] = 2 * S[0] - 2.0 / (S[0] * S[0] * S[0]);
            es[1] = 2 * S[1] - 2.0 / (S[1] * S[1] * S[1]);
            break;
            
        case 3: // amips
            es[0] = amips_param * exp(amips_param * (S[0] / S[1] + S[1] / S[0])) * (1.0 / S[1] - S[1] / (S[0] * S[0]));
            es[1] = amips_param * exp(amips_param * (S[0] / S[1] + S[1] / S[0])) * (1.0 / S[0] - S[0] / (S[1] * S[1]));
            break;
            
        case 4: // conf
            es[0] = 2 * S[0] / (S[1] * S[1]);
            es[1] = -2 * S[0] * S[0] / (S[1] * S[1] * S[1]);
            break;
            
        case 5: // gmr
            es[0] = c1_param * (1.0 / S[1] - S[1] / (S[0] * S[0])) - d1_param *2.*S[1]*(1.0 - S[0]*S[1]);
            es[1] = c1_param * (1.0 / S[0] - S[0] / (S[1] * S[1])) - d1_param *2.*S[0]*(1.0 - S[0]*S[1]);

            break;
				
        case 6: // olg
            
            if(S[0] > 1.0 / S[1])
            {
                es[0] = 1;
                es[1] = 0;
            }else{
                es[0] = 0;
                es[1] = -1.0 / (S[1] * S[1]);
            }
            
            break;
    }
    
    return es;

}


double energy_hessian(double s1, double s2, double s1j, double s1k, double s2j, double s2k, double s1jk, double s2jk, int type)
{
    
    double hessian = 0;
    double v1, v2;
    
    switch(type)
    {
        case 0: //arap
            v1 = 2 * ((s1 - 1) * s1jk + (s2 - 1) * s2jk);
            v2 = 2 * (s1j * s1k + s2j * s2k);
            hessian = v1 + v2;
            break;
            
        case 1: // mips
            
            v1 = - s2k / pow(s2, 2) - s2k /pow(s1, 2) + 2 * s2 * s1k / pow(s1, 3);
            v2 = - s1k / pow(s1, 2) - s1k /pow(s2, 2) + 2 * s1 * s2k / pow(s2, 3);
            
            hessian = v1 * s1j + v2 * s2j;
            
            v1 = 1.0 / s2 - s2 / pow(s1, 2);
            v2 = 1.0 / s1 - s1 / pow(s2, 2);
            
            hessian += v1 * s1jk + v2 * s2jk;
            
            break;
            
        case 2: // iso
            
            v1 = 2 * s1k * s1j + 2 * s1 * s1jk;
            v2 = 2 * s2k * s2j + 2 * s2 * s2jk;
            
            v1 += 6 * s1k * s1j / pow(s1, 4) - 2 * s1jk / pow(s1, 3);
            v2 += 6 * s2k * s2j / pow(s2, 4) - 2 * s2jk / pow(s2, 3);
            
            hessian = v1 + v2;
            
            break;
            
        case 3: // amips
            
            break;
            
        case 4: // conf
            v1 = (2 / (s1 * s1) * s1k - 4 * s1 / (s2 * s2 * s2) * s2k) * s1j + 2 * s1 / (s2 * s2) * s1jk;
            v2 = (6 * s1 * s1 / (s2 * s2 * s2 * s2) * s2k - 4 * s1 / (s2 * s2 * s2) * s1k) * s2j - 2 * s1 * s1 / (s2 * s2 * s2) * s2jk;
            
            hessian = v1 + v2;
            break;
            
        case 5:
            
            double e1 = c1_param * (1 / s2 - s2 / (s1 * s1)) + 2 * d1_param * s2 * (s1 * s2 - 1);
            double e2 = c1_param * (1 / s1 - s1 / (s2 * s2)) + 2 * d1_param * s1 * (s1 * s2 - 1);
            double e11 = 2 * c1_param * s2 / (s1 * s1 * s1) + 2 * d1_param * s2 * s2;
            double e22 = 2 * c1_param * s1 / (s2 * s2 * s2) + 2 * d1_param * s1 * s1;
            double e12 = c1_param * (-1 / (s1 * s1) - 1 / (s2 * s2)) + 2 * d1_param * (2 * s1 * s2 - 1);
            
            hessian = (e11 * s1k + e12 * s2k) * s1j + e1 * s1jk + (e12 * s1k + e22 * s2k) * s2j + e2 * s2jk;
            
            break;
    }
    
    return hessian;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *row_mex, *col_mex, *val_mex, *grad_mex, *J_value_mex, *JT_value_mex, *wu_mex, *bu_mex;
    
    double *tri_num, *X_g_inv, *tri_areas, *obj_tri, *q_target, *type, *amips_s, *F_dot, *ver_num, *c1_g, *c2_g, *d1_g, *clamp, *J_index, *JT_index, *JTJ_info, *x2u;
    double *row, *col, *val, *grad, *J_value, *JT_value, *wu, *bu;
    

    // input list
    tri_num = mxGetPr(prhs[0]);
    X_g_inv = mxGetPr(prhs[1]);
    tri_areas = mxGetPr(prhs[2]);
    obj_tri = mxGetPr(prhs[3]);
    q_target = mxGetPr(prhs[4]);
    type = mxGetPr(prhs[5]);
    amips_s = mxGetPr(prhs[6]);
    F_dot = mxGetPr(prhs[7]);
    ver_num = mxGetPr(prhs[8]);
    c1_g = mxGetPr(prhs[9]);
    c2_g = mxGetPr(prhs[10]);
    d1_g = mxGetPr(prhs[11]);
    clamp = mxGetPr(prhs[12]);
    J_index = mxGetPr(prhs[13]);
    JT_index = mxGetPr(prhs[14]);
    JTJ_info = mxGetPr(prhs[15]);
    x2u = mxGetPr(prhs[16]);
    

    int tri_n = tri_num[0];
    int ver_n = ver_num[0];
    int energy_type = type[0];
    double proj_threshold = clamp[0];
    amips_param = amips_s[0];
    c1_param = c1_g[0];
    c2_param = c2_g[0];
    d1_param = d1_g[0];
    int J_rn = JTJ_info[0];
    int J_cn = JTJ_info[1];
    int JT_rn = JTJ_info[2];
    int JT_cn = JTJ_info[3];
    

    // output list
    grad_mex = plhs[0] = mxCreateDoubleMatrix(2 * ver_n, 1, mxREAL);
    J_value_mex = plhs[1] = mxCreateDoubleMatrix(J_rn * J_cn, 1, mxREAL);
    JT_value_mex = plhs[2] = mxCreateDoubleMatrix(JT_rn * JT_cn, 1, mxREAL);
    wu_mex = plhs[3] = mxCreateDoubleMatrix(tri_n, 1, mxREAL);
    bu_mex = plhs[4] = mxCreateDoubleMatrix(tri_n, 1, mxREAL);
    row_mex = plhs[5] = mxCreateDoubleMatrix(36 * tri_n, 1, mxREAL);
    col_mex = plhs[6] = mxCreateDoubleMatrix(36 * tri_n, 1, mxREAL);
    val_mex = plhs[7] = mxCreateDoubleMatrix(36 * tri_n, 1, mxREAL);
    
    
    grad = mxGetPr(grad_mex);
    J_value = mxGetPr(J_value_mex);
    JT_value = mxGetPr(JT_value_mex);
    wu = mxGetPr(wu_mex);
    bu = mxGetPr(bu_mex);
    row = mxGetPr(row_mex);
    col = mxGetPr(col_mex);
    val = mxGetPr(val_mex);
    
    memset(J_value, 0, J_rn * J_cn * sizeof(double));
    memset(JT_value, 0, JT_rn * JT_cn * sizeof(double));
    memset(grad, 0, 2 * ver_n * sizeof(double));

    vector<int>* offset = new vector<int>;
    
    offset->resize(2 * ver_n, 0);
    
    int tri[3];
    double tri_area;
    int index[6];
    
    Vector2d es;
    Vector2d cs;
    
    double es1, es2;
    double d1, d2;
    double dphi;
    
    Matrix2d B, X_f, A, F_d, T_ddot, T;
    Vector2d S, u0, u1, v0, v1, temp_w;
    Matrix2d U, V;
    MatrixXd ELS(6, 6);
    VectorXd ELS_S;
    MatrixXd ELS_X;
    DiagonalMatrix<double, 6> ELS_D;
    
    Matrix2d T_dot[6];
    Matrix2d w_U_dot[6];
    Matrix2d w_V_dot[6];
    
    Matrix2d temp_A, inv_A;
    
    Matrix2d II;
    II(0, 0) = 1;
    II(1, 1) = 1;
    II(0, 1) = 0;
    II(1, 0) = 0;
    
    double tao = 0.01;
    
    double value;
 
    for(int i = 0; i < tri_n; i++)
    {
        //mexPrintf("%d %d %d\n", (int)obj_tri[i], (int)obj_tri[i + tri_n], (int)obj_tri[i + 2 * tri_n]);
        tri[0] = obj_tri[i] - 1; 
        tri[1] = obj_tri[i + tri_n] - 1; 
        tri[2] = obj_tri[i + 2 * tri_n] - 1; 
        index[0] = 2 * tri[0];
        index[1] = index[0] + 1;
        index[2] = 2 * tri[1];
        index[3] = index[2] + 1;
        index[4] = 2 * tri[2];
        index[5] = index[4] + 1;
        tri_area = tri_areas[i];
        
        B(0, 0) = X_g_inv[i];
        B(0, 1) = X_g_inv[i + tri_n * 2];
        B(1, 0) = X_g_inv[i + tri_n];
        B(1, 1) = X_g_inv[i + tri_n * 3];
        
        X_f(0, 0) = q_target[2 * tri[1]] - q_target[2 * tri[0]];
        X_f(1, 0) = q_target[2 * tri[1] + 1] - q_target[2 * tri[0] + 1];  
        X_f(0, 1) = q_target[2 * tri[2]] - q_target[2 * tri[0]];
        X_f(1, 1) = q_target[2 * tri[2] + 1] - q_target[2 * tri[0] + 1];
        
        A = X_f * B;
        
        JacobiSVD<Matrix2d> svd(A, ComputeFullU | ComputeFullV);
        
        S = svd.singularValues();
        U = svd.matrixU();
        V = svd.matrixV();
        
        bu[i] = S[0] * S[1];
        
        es = energy_derivative(energy_type, S);
        
        cs = S;
        
        es1 = es[0];
        es2 = es[1];
        
        /////////////////////////////////////////////////////////////////////////////////////
        
        T(0, 0) = S[0];
        T(0, 1) = 0;
        T(1, 0) = 0;
        T(1, 1) = S[1];
        
        temp_A(0, 0) = S[1];
        temp_A(0, 1) = S[0];
        temp_A(1, 0) = S[0];
        temp_A(1, 1) = S[1];
        
        u0 = U.col(0);
        u1 = U.col(1);
        v0 = V.col(0);
        v1 = V.col(1);
        
        if(abs(S[0] - S[1]) < 1e-5)
        {
            inv_A = temp_A.transpose() * temp_A + tao * II;
            inv_A = inv_A.inverse() * temp_A.transpose();
            //mexPrintf("ill conditioned\n");
        }else{
            inv_A = temp_A.inverse();
        }
        
        wu[i] = 0;
        
        for(int j = 0; j < 6; j++)
        {
            F_d(0, 0) = F_dot[i + j * tri_n];
            F_d(0, 1) = F_dot[i + j * tri_n + 6 * 2 * tri_n];
            F_d(1, 0) = F_dot[i + j * tri_n + 6 * tri_n];
            F_d(1, 1) = F_dot[i + j * tri_n + 6 * 3 * tri_n];
            
            d1 = u0.transpose() * F_d * v0;
            d2 = u1.transpose() * F_d * v1;
            
            dphi = d1 * es1 + d2 * es2;
            
            int cache_1 = 2 * tri[j / 2] + (j % 2);
            int cache_2 = x2u[cache_1];
            
            grad[cache_1] += tri_area * dphi;
            
            if(cache_2 > 0)
            {
                int idx = cache_2 - 1;
                
                double value = d1 * cs[1] + d2 * cs[0];
                
                J_value[idx * J_cn + (*offset)[cache_1]] = value;
                
                (*offset)[cache_1]++;
                
                JT_value[i * JT_cn + j] = value;
                
                wu[i] += value * value;
            }
               
            T_dot[j](0, 0) = d1;
            T_dot[j](0, 1) = 0;
            T_dot[j](1, 0) = 0;
            T_dot[j](1, 1) = d2;
            
            temp_w[0] = u0.transpose() * F_d * v1;
            temp_w[1] = -1.0 * u1.transpose() * F_d * v0;
            temp_w = inv_A * temp_w;
            
            w_U_dot[j](0, 0) = 0;
            w_U_dot[j](0, 1) = temp_w[0];
            w_U_dot[j](1, 0) = -1 * temp_w[0];
            w_U_dot[j](1, 1) = 0;
            
            w_V_dot[j](0, 0) = 0;
            w_V_dot[j](0, 1) = temp_w[1];
            w_V_dot[j](1, 0) = -1 * temp_w[1];
            w_V_dot[j](1, 1) = 0;           
        }
        
        for(int j = 0; j < 6; j++)
        {
            for(int k = j; k < 6; k++)
            {
                T_ddot = w_U_dot[k].transpose() * w_U_dot[j] * T + T * w_V_dot[j] * w_V_dot[k].transpose() - w_U_dot[j] * T * w_V_dot[k] - w_U_dot[k] * T * w_V_dot[j];
                
                value = energy_hessian(S[0], S[1], T_dot[j](0, 0), T_dot[k](0, 0), T_dot[j](1, 1), T_dot[k](1, 1), T_ddot(0, 0), T_ddot(1, 1), energy_type);
                
                ELS(j, k) = tri_area * value;
                ELS(k, j) = tri_area * value;
                
                row[36 * i + j * 6 + k] = index[j];
                col[36 * i + j * 6 + k] = index[k];
                
                row[36 * i + k * 6 + j] = index[k];
                col[36 * i + k * 6 + j] = index[j];
            }
        }
        
        SelfAdjointEigenSolver<MatrixXd> eig(ELS);
        
        ELS_S = eig.eigenvalues();
        
        ELS_X = eig.eigenvectors();
        
        for(int j = 0; j < 6; j++)
        {
            ELS_S[j] = fmax(ELS_S[j], proj_threshold);
        }
        
        ELS_D.diagonal() << ELS_S[0], ELS_S[1], ELS_S[2], ELS_S[3], ELS_S[4], ELS_S[5];
        
        ELS = ELS_X * ELS_D * ELS_X.transpose();
        
        for(int j = 0; j < 6; j++)
        {
            for(int k = j; k < 6; k++)
            {         
                val[36 * i + j * 6 + k] = ELS(j, k);
                val[36 * i + k * 6 + j] = ELS(k, j);
            }
        }
      
    }
    
    offset->clear();
    delete offset;
    
    return;
    
}