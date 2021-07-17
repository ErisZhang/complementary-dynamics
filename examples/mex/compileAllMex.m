eigenFolder = '../../libigl/external/eigen';
iglFolder = '../../libigl/include';

mex(['-I',eigenFolder],'energy_hessian_mex.cpp');
mex(['-I',eigenFolder],'precompute_mex.cpp');
mex(['-I',eigenFolder],'energy_value_mex.cpp');
mex(['-I',eigenFolder],['-I',iglFolder],'read_3d_bone_anim.cpp');
mex(['-I',eigenFolder],['-I',iglFolder],'read_3d_pnt_anim.cpp');
mex(['-I',eigenFolder],['-I',iglFolder],'read_blendshape_anim.cpp');
