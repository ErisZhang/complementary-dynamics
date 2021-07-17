#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>


using namespace Eigen;
using namespace std;


double amips_param;
double c1_param;
double c2_param;
double d1_param;


double energy_value(int type, Vector2d S)
{
    
    double value = 0;
    double J;
    double l1;
    double l2;
    
    switch(type)
    {
        case 0: //arap
            value = (S[0] - 1) * (S[0] - 1) + (S[1] - 1) * (S[1] - 1);
            break;
            
        case 1: // mips
            value = S[0] / S[1] + S[1] / S[0];
            //value *= 100;
            break;
            
        case 2: // iso
            value = S[0] * S[0] + 1.0 / (S[0] * S[0]) + S[1] * S[1] + 1.0 / (S[1] * S[1]);
            break;
            
        case 3: // amips
            value = exp(amips_param * (S[0] / S[1] + S[1] / S[0]));
            break;
            
        case 4: // conf
            value = S[0] / S[1];
            value *= value;
            break;
            
        case 5:

						//            J = S[0] * S[1];
						//            l1 = S[0] * S[0] + S[1] * S[1];
						//            l2 = S[0] * S[0] * S[1] * S[1];
						//            //value = (J - 1) * (J - 1);
						//            value = c1_param * (pow(J, -2.0 / 3.0) * l1 - 3) + c2_param * (pow(J, -4.0 / 3.0) * l2 - 3) + d1_param * (J - 1) * (J - 1);
						//            break;
				
						J = S[0] * S[1];
						value = c1_param * (S[0] / S[1] + S[1] / S[0] - 2) + d1_param * (J - 1) * (J - 1);
						break;
				
        case 6: // olg
				
            value = max(S[0], 1.0 / S[1]);
            break;
    }
    
    return value;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex;
    const int *dims;
    double *tri_num, *X_g_inv, *tri_areas, *obj_tri, *q_target, *type, *amips_s, *c1, *c2, *d1;
    double *output;
    
    tri_num = mxGetPr(prhs[0]);
    X_g_inv = mxGetPr(prhs[1]);
    tri_areas = mxGetPr(prhs[2]);
    obj_tri = mxGetPr(prhs[3]);
    q_target = mxGetPr(prhs[4]);
    type = mxGetPr(prhs[5]);
    amips_s = mxGetPr(prhs[6]);
    c1 = mxGetPr(prhs[7]);
    c2 = mxGetPr(prhs[8]);
    d1 = mxGetPr(prhs[9]);
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    output = mxGetPr(output_mex);
    
    output[0] = 0;
    
    int tri_n = tri_num[0];
    int energy_type = type[0];
    amips_param = amips_s[0];
    c1_param = c1[0];
    c2_param = c2[0];
    d1_param = d1[0];
    
    int tri[3];
    double tri_area;
    
    Matrix2d B, X_f, A;
    Vector2d S;
    Matrix2d U, V;
 
    for(int i = 0; i < tri_n; i++)
    {
        //mexPrintf("%d %d %d\n", (int)obj_tri[i], (int)obj_tri[i + tri_n], (int)obj_tri[i + 2 * tri_n]);
        tri[0] = obj_tri[i] - 1; 
        tri[1] = obj_tri[i + tri_n] - 1; 
        tri[2] = obj_tri[i + 2 * tri_n] - 1; 
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
        
        double cache = A.determinant();
        
        if(cache <= 0)
        {
            S[1] *= -1;
        }
        
        double deviation = energy_value(energy_type, S);
        
        if(cache <= 0 && energy_type != 0)
        {
            mexPrintf("inverted\n");
            output[0] = 1e20;
            break;
        }
        
        output[0] += tri_area * deviation;
        
    }
    
    return;
    
}