function R = transform2d(theta,x,y)
    R=zeros(3,3);
    cost = cos(theta);
    sint = sin(theta);
    R(1,1) = cost;
    R(1,2) = -sint;
    R(2,1) = sint;
    R(2,2) = cost;
    R(1,3) = -x*cost+y*sint+x;
    R(2,3) = -x*sint-y*cost+y;
    R(3,3) = 1;
end