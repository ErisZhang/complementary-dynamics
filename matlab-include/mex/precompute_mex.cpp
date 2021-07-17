#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>


using namespace Eigen;
using namespace std;


double get_tri_area(Vector2d p0, Vector2d p1, Vector2d p2)
{
    
    Vector2d u = p1 - p0;
    Vector2d v = p2 - p0;
    
    return 0.5 * sqrt(u.dot(u) * v.dot(v) - u.dot(v) * u.dot(v));
    
}


double get_tri_area(Vector3d p0, Vector3d p1, Vector3d p2)
{
    
    Vector3d u = p1 - p0;
    Vector3d v = p2 - p0;
    
    return 0.5 * sqrt(u.dot(u) * v.dot(v) - u.dot(v) * u.dot(v));
    
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *X_g_inv_mex, *tri_areas_mex, *F_dot_mex, *row_mex, *col_mex, *val_mex, *x2u_mex, *J_mex, *J_info_mex, *JT_mex, *JT_info_mex, *perimeter_mex;
    
    double *tri_num, *tri_list, *ver_num, *ver_list, *uv_mesh, *dirichlet;
    double *X_g_inv, *tri_areas, *F_dot, *row, *col, *val, *x2u, *J, *J_info, *JT, *JT_info, *perimeter;
    
    tri_num = mxGetPr(prhs[0]);
    tri_list = mxGetPr(prhs[1]);
    ver_num = mxGetPr(prhs[2]);
    ver_list = mxGetPr(prhs[3]);
    uv_mesh = mxGetPr(prhs[4]);
    dirichlet = mxGetPr(prhs[5]);
    
    int tri_n = tri_num[0];
    int ver_n = ver_num[0];
    int is_uv_mesh = uv_mesh[0];
    
    X_g_inv_mex = plhs[0] = mxCreateDoubleMatrix(4 * tri_n, 1, mxREAL);
    tri_areas_mex = plhs[1] = mxCreateDoubleMatrix(tri_n, 1, mxREAL);
    F_dot_mex = plhs[2] = mxCreateDoubleMatrix(24 * tri_n, 1, mxREAL);
    row_mex = plhs[3] = mxCreateDoubleMatrix(18 * tri_n, 1, mxREAL);
    col_mex = plhs[4] = mxCreateDoubleMatrix(18 * tri_n, 1, mxREAL);
    val_mex = plhs[5] = mxCreateDoubleMatrix(18 * tri_n, 1, mxREAL);
    x2u_mex = plhs[6] = mxCreateDoubleMatrix(2 * ver_n, 1, mxREAL);
    
    X_g_inv = mxGetPr(X_g_inv_mex);
    tri_areas = mxGetPr(tri_areas_mex);
    F_dot = mxGetPr(F_dot_mex);
    row = mxGetPr(row_mex);
    col = mxGetPr(col_mex);
    val = mxGetPr(val_mex);
    x2u = mxGetPr(x2u_mex);
    
    int count = 1;
    int fixed_num = 0;
    int tri[3];
    
    for(int i = 0; i < 2 * ver_n; i++)
    {        
        if(dirichlet[i] == 0)
        {
            x2u[i] = count;
            count++;
        }else
        {
            fixed_num++;
            x2u[i] = 0;
        }
    }
    
    //mexPrintf("check 1\n");
    
    vector<vector<int>*> ver_tri_map; 
    
    for(int i = 0; i < ver_n; i++)
    {
        vector<int>* tris = new vector<int>();
        ver_tri_map.push_back(tris);
    }
    
    for(int i = 0; i < tri_n; i++)
    {
        tri[0] = tri_list[i] - 1; 
        tri[1] = tri_list[i + tri_n] - 1; 
        tri[2] = tri_list[i + 2 * tri_n] - 1;
        
        ver_tri_map[tri[0]]->push_back(i);
        ver_tri_map[tri[1]]->push_back(i);
        ver_tri_map[tri[2]]->push_back(i);
    }
    
    int max_valence = 0;
    
    for(int i = 0; i < ver_n; i++)
    {
        if(max_valence < ver_tri_map[i]->size())
        {
            max_valence = ver_tri_map[i]->size();
        }
    }
    
    J_mex = plhs[7] = mxCreateDoubleMatrix((2 * ver_n - fixed_num) * max_valence, 1, mxREAL);
    J_info_mex = plhs[8] = mxCreateDoubleMatrix(2, 1, mxREAL);
    JT_mex = plhs[9] = mxCreateDoubleMatrix(tri_n * 3 * 2, 1, mxREAL);
    JT_info_mex = plhs[10] = mxCreateDoubleMatrix(2, 1, mxREAL);
    perimeter_mex = plhs[11] = mxCreateDoubleMatrix(ver_n, 1, mxREAL);
    
    J = mxGetPr(J_mex);
    J_info = mxGetPr(J_info_mex);
    JT = mxGetPr(JT_mex);
    JT_info = mxGetPr(JT_info_mex);
    perimeter = mxGetPr(perimeter_mex);
    
    count = 0;
    
    J_info[0] = 2 * ver_n - fixed_num;
    J_info[1] = max_valence;
    
    for(int i = 0; i < ver_n; i++)
    {
        perimeter[i] = 0;
        
        if(dirichlet[2 * i] == 0)
        {
            for(int j = 0; j < max_valence; j++)
            {
                if(j < ver_tri_map[i]->size())
                {
                    J[count * max_valence + j] = ver_tri_map[i]->at(j);
                }else{
                    J[count * max_valence + j] = -1;
                }
            }
            
            count++;
            
            for(int j = 0; j < max_valence; j++)
            {
                if(j < ver_tri_map[i]->size())
                {
                    J[count * max_valence + j] = ver_tri_map[i]->at(j);
                }else{
                    J[count * max_valence + j] = -1;
                }
            }
            
            count++;
        }
    }
    
    
    JT_info[0] = tri_n;
    JT_info[1] = 3 * 2;
    
    for(int i = 0; i < tri_n; i++)
    { 
        tri[0] = tri_list[i] - 1; 
        tri[1] = tri_list[i + tri_n] - 1; 
        tri[2] = tri_list[i + 2 * tri_n] - 1;
        
        for(int j = 0; j < 3; j++)
        {
            if(x2u[2 * tri[j]] > 0)
            {
                JT[i * 2 * 3 + 2 * j + 0] = x2u[2 * tri[j] + 0] - 1;
                JT[i * 2 * 3 + 2 * j + 1] = x2u[2 * tri[j] + 1] - 1;
            }else{
                JT[i * 2 * 3 + 2 * j + 0] = -1;
                JT[i * 2 * 3 + 2 * j + 1] = -1;
            }     
        }
    }
    
    for(int i = 0; i < ver_n; i++)
    {
        ver_tri_map[i]->clear();
        delete ver_tri_map[i];
    }
    
    ver_tri_map.clear();
    
    ///////////////////////////////////////////////////////////////////////
    
    Matrix2d X_g, B;
    
    double C_sub[4][6];
    
    double tri_area;
        
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            C_sub[i][j] = 0;
        }
    }
  
    for(int i = 0; i < tri_n; i++)
    {
        
        tri[0] = tri_list[i] - 1; 
        tri[1] = tri_list[i + tri_n] - 1; 
        tri[2] = tri_list[i + 2 * tri_n] - 1; 
        
        /////////////////////////////////////////////////////////////
        
        if(is_uv_mesh)
        {
            
            Vector3d p0(ver_list[3 * tri[0]], ver_list[3 * tri[0] + 1], ver_list[3 * tri[0] + 2]);
            Vector3d p1(ver_list[3 * tri[1]], ver_list[3 * tri[1] + 1], ver_list[3 * tri[1] + 2]);
            Vector3d p2(ver_list[3 * tri[2]], ver_list[3 * tri[2] + 1], ver_list[3 * tri[2] + 2]);
            
            perimeter[tri[0]] += (p1 - p2).norm();
            perimeter[tri[1]] += (p2 - p0).norm();
            perimeter[tri[2]] += (p0 - p1).norm();
            
            tri_area = tri_areas[i] = get_tri_area(p0, p1, p2);
            
            Vector3d bx = p1 - p0;
            Vector3d cx = p2 - p0;
            
            Vector3d Ux = bx;
            Ux.normalize();
            
            Vector3d w = Ux.cross(cx);
            Vector3d Wx = w;
            Wx.normalize();
            
            Vector3d Vx = Ux.cross(Wx);
            
            Matrix3d R;
            R.col(0) = Ux;
            R.col(1) = Vx;
            R.col(2) = Wx;
            
            Vector3d vb = R.transpose() * bx;
            Vector3d vc = R.transpose() * cx;
            
            X_g << vb[0], vc[0], vb[1], vc[1];
            
            if(X_g.determinant() < 0)
            {
                X_g.row(1) = -1.0 * X_g.row(1);
            }
            
        }else{
            
            Vector2d p0(ver_list[2 * tri[0]], ver_list[2 * tri[0] + 1]);
            Vector2d p1(ver_list[2 * tri[1]], ver_list[2 * tri[1] + 1]);
            Vector2d p2(ver_list[2 * tri[2]], ver_list[2 * tri[2] + 1]);
            
            perimeter[tri[0]] += (p1 - p2).norm();
            perimeter[tri[1]] += (p2 - p0).norm();
            perimeter[tri[2]] += (p0 - p1).norm();
        
            tri_area = tri_areas[i] = get_tri_area(p0, p1, p2);
        
            Vector2d e0 = p1 - p0;
            Vector2d e1 = p2 - p0;
        
            X_g << e0[0], e1[0], e0[1], e1[1];
            
        }
          
        ///////////////////////////////////////////////////////////////
        
        B = X_g.inverse();
        
        X_g_inv[i] = B(0, 0);
        X_g_inv[tri_n + i] = B(1, 0);
        X_g_inv[2 * tri_n + i] = B(0, 1);
        X_g_inv[3 * tri_n + i] = B(1, 1);
        
        C_sub[0][0] = -(B(0, 0) + B(1, 0));
        C_sub[1][1] = -(B(0, 0) + B(1, 0));
        C_sub[0][2] = B(0, 0);
        C_sub[1][3] = B(0, 0);
        C_sub[0][4] = B(1, 0);
        C_sub[1][5] = B(1, 0);
            
        C_sub[2][0] = -(B(0, 1) + B(1, 1));
        C_sub[3][1] = -(B(0, 1) + B(1, 1));
        C_sub[2][2] = B(0, 1); 
        C_sub[3][3] = B(0, 1);
        C_sub[2][4] = B(1, 1);
        C_sub[3][5] = B(1, 1);
        
        for(int j = 0; j < 6; j++)
        {
            F_dot[0 * tri_n * 6 + j * tri_n + i] = C_sub[0][j];
            F_dot[1 * tri_n * 6 + j * tri_n + i] = C_sub[1][j];
            F_dot[2 * tri_n * 6 + j * tri_n + i] = C_sub[2][j];
            F_dot[3 * tri_n * 6 + j * tri_n + i] = C_sub[3][j];
        }
        
        //////////////////////////////////////////////////////////////////////
        
        Matrix2d A0;
        A0 << -1.0 * (B(0, 0) + B(1, 0)), 0, 0, -1.0 * (B(0, 0) + B(1, 0));
        Matrix2d A1;
        A1 << -1.0 * (B(0, 1) + B(1, 1)), 0, 0, -1.0 * (B(0, 1) + B(1, 1));
        Matrix2d B0;
        B0 << B(0, 0), 0, 0, B(0, 0);
        Matrix2d B1;
        B1 << B(0, 1), 0, 0, B(0, 1);
        Matrix2d C0;
        C0 << B(1, 0), 0, 0, B(1, 0);
        Matrix2d C1;
        C1 << B(1, 1), 0, 0, B(1, 1);
        
        Matrix2d AA = A0 * A0 + A1 * A1;
        Matrix2d AB = A0 * B0 + A1 * B1;
        Matrix2d AC = A0 * C0 + A1 * C1;
        Matrix2d BA = B0 * A0 + B1 * A1;
        Matrix2d BB = B0 * B0 + B1 * B1;
        Matrix2d BC = B0 * C0 + B1 * C1;
        Matrix2d CA = C0 * A0 + C1 * A1;
        Matrix2d CB = C0 * B0 + C1 * B1;
        Matrix2d CC = C0 * C0 + C1 * C1;

        row[18 * i + 0] = 2 * tri[0];
        col[18 * i + 0] = 2 * tri[0];
        val[18 * i + 0] = AA(0, 0) * tri_area;
        
        row[18 * i + 1] = 2 * tri[0] + 1;
        col[18 * i + 1] = 2 * tri[0] + 1;
        val[18 * i + 1] = AA(1, 1) * tri_area;
        
        row[18 * i + 2] = 2 * tri[0];
        col[18 * i + 2] = 2 * tri[1];
        val[18 * i + 2] = AB(0, 0) * tri_area;
        
        row[18 * i + 3] = 2 * tri[0] + 1;
        col[18 * i + 3] = 2 * tri[1] + 1;
        val[18 * i + 3] = AB(1, 1) * tri_area;
        
        row[18 * i + 4] = 2 * tri[0];
        col[18 * i + 4] = 2 * tri[2];
        val[18 * i + 4] = AC(0, 0) * tri_area;
        
        row[18 * i + 5] = 2 * tri[0] + 1;
        col[18 * i + 5] = 2 * tri[2] + 1;
        val[18 * i + 5] = AC(1, 1) * tri_area;
        
        row[18 * i + 6] = 2 * tri[1];
        col[18 * i + 6] = 2 * tri[0];
        val[18 * i + 6] = BA(0, 0) * tri_area;
        
        row[18 * i + 7] = 2 * tri[1] + 1;
        col[18 * i + 7] = 2 * tri[0] + 1;
        val[18 * i + 7] = BA(1, 1) * tri_area;
        
        row[18 * i + 8] = 2 * tri[1];
        col[18 * i + 8] = 2 * tri[1];
        val[18 * i + 8] = BB(0, 0) * tri_area;
        
        row[18 * i + 9] = 2 * tri[1] + 1;
        col[18 * i + 9] = 2 * tri[1] + 1;
        val[18 * i + 9] = BB(1, 1) * tri_area;
        
        row[18 * i + 10] = 2 * tri[1];
        col[18 * i + 10] = 2 * tri[2];
        val[18 * i + 10] = BC(0, 0) * tri_area;
        
        row[18 * i + 11] = 2 * tri[1] + 1;
        col[18 * i + 11] = 2 * tri[2] + 1;
        val[18 * i + 11] = BC(1, 1) * tri_area;
        
        row[18 * i + 12] = 2 * tri[2];
        col[18 * i + 12] = 2 * tri[0];
        val[18 * i + 12] = CA(0, 0) * tri_area;
        
        row[18 * i + 13] = 2 * tri[2] + 1;
        col[18 * i + 13] = 2 * tri[0] + 1;
        val[18 * i + 13] = CA(1, 1) * tri_area;
        
        row[18 * i + 14] = 2 * tri[2];
        col[18 * i + 14] = 2 * tri[1];
        val[18 * i + 14] = CB(0, 0) * tri_area;
        
        row[18 * i + 15] = 2 * tri[2] + 1;
        col[18 * i + 15] = 2 * tri[1] + 1;
        val[18 * i + 15] = CB(1, 1) * tri_area;
        
        row[18 * i + 16] = 2 * tri[2];
        col[18 * i + 16] = 2 * tri[2];
        val[18 * i + 16] = CC(0, 0) * tri_area;
        
        row[18 * i + 17] = 2 * tri[2] + 1;
        col[18 * i + 17] = 2 * tri[2] + 1;
        val[18 * i + 17] = CC(1, 1) * tri_area;
        
    }
    
    return;
    
}