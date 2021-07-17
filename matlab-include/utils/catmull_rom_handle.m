function [X,Y,tX,tY,nX,nY] = catmull_rom_handle(handleX,handleY,num)
    
    assert(size(handleX,2) > 1 || size(handleX,1) > 1);
    handleX = reshape(handleX,1,size(handleX,1)*size(handleX,2));
    handleY = reshape(handleY,1,size(handleY,1)*size(handleY,2));

    T = linspace(0,1,(size(handleX,2)-1)*num+1);
    handleT = T(1:num:end);
    X = zeros(size(T));
    Y = zeros(size(T));
    tX = zeros(size(T));
    tY = zeros(size(T));
    nX = zeros(size(T));
    nY = zeros(size(T));

    for i = 1:(size(T,2)-1)
        idx0 = ceil(i/num);
        idx1 = ceil(i/num)+1;
        t0 = handleT(idx0);
        t1 = handleT(idx1);
        t = T(i);
        t = (t-t0)/(t1-t0);
        if idx0 == 1 && idx1 == size(handleT,2)
            mX0 = (handleX(idx0+1)-handleX(idx0))/(handleT(idx0+1)-handleT(idx0));
            mX1 = (handleX(idx1)-handleX(idx1-1))/(handleT(idx1)-handleT(idx1-1));
            mY0 = (handleY(idx0+1)-handleY(idx0))/(handleT(idx0+1)-handleT(idx0));
            mY1 = (handleY(idx1)-handleY(idx1-1))/(handleT(idx1)-handleT(idx1-1));
        elseif idx0 == 1 && idx1 ~= size(handleT,2)
            mX0 = (handleX(idx0+1)-handleX(idx0))/(handleT(idx0+1)-handleT(idx0));
            mX1 = (handleX(idx1+1)-handleX(idx1-1))/(handleT(idx1+1)-handleT(idx1-1));
            mY0 = (handleY(idx0+1)-handleY(idx0))/(handleT(idx0+1)-handleT(idx0));
            mY1 = (handleY(idx1+1)-handleY(idx1-1))/(handleT(idx1+1)-handleT(idx1-1));
        elseif idx0 ~= 1 && idx1 == size(handleT,2)
            mX0 = (handleX(idx0+1)-handleX(idx0-1))/(handleT(idx0+1)-handleT(idx0-1));
            mX1 = (handleX(idx1)-handleX(idx1-1))/(handleT(idx1)-handleT(idx1-1));
            mY0 = (handleY(idx0+1)-handleY(idx0-1))/(handleT(idx0+1)-handleT(idx0-1));
            mY1 = (handleY(idx1)-handleY(idx1-1))/(handleT(idx1)-handleT(idx1-1));
        else
            mX0 = (handleX(idx0+1)-handleX(idx0-1))/(handleT(idx0+1)-handleT(idx0-1));
            mX1 = (handleX(idx1+1)-handleX(idx1-1))/(handleT(idx1+1)-handleT(idx1-1));
            mY0 = (handleY(idx0+1)-handleY(idx0-1))/(handleT(idx0+1)-handleT(idx0-1));
            mY1 = (handleY(idx1+1)-handleY(idx1-1))/(handleT(idx1+1)-handleT(idx1-1));
        end
        pX0 = handleX(idx0);
        pX1 = handleX(idx1);
        pY0 = handleY(idx0);
        pY1 = handleY(idx1);
        pX = (2*t^3-3*t^2+1)*pX0 + (t^3 - 2*t^2 + t)*(handleT(idx1)-handleT(idx0))*mX0 + (-2*t^3+3*t^2)*pX1 + (t^3 -t^2)*(handleT(idx1)-handleT(idx0))*mX1;
        pY = (2*t^3-3*t^2+1)*pY0 + (t^3 - 2*t^2 + t)*(handleT(idx1)-handleT(idx0))*mY0 + (-2*t^3+3*t^2)*pY1 + (t^3 -t^2)*(handleT(idx1)-handleT(idx0))*mY1;
        ptX = (6*t^2-6*t)*pX0 + (3*t^2 - 4*t + 1)*(handleT(idx1)-handleT(idx0))*mX0 + (-6*t^2+6*t)*pX1 + (3*t^2 -2*t)*(handleT(idx1)-handleT(idx0))*mX1;
        ptY = (6*t^2-6*t)*pY0 + (3*t^2 - 4*t + 1)*(handleT(idx1)-handleT(idx0))*mY0 + (-6*t^2+6*t)*pY1 + (3*t^2 -2*t)*(handleT(idx1)-handleT(idx0))*mY1;
        N = [0 -1; 1 0]*[ptX;ptY];
        X(i) = pX;
        Y(i) = pY;
        tX(i) = ptX/norm([ptX ptY]);
        tY(i) = ptY/norm([ptX ptY]);
        nX(i) = N(1)/norm(N);
        nY(i) = N(2)/norm(N);
    end
    
    % last handle
    X(end) = handleX(end);
    Y(end) = handleY(end);
    tX(end) = mX1/norm([mX1 mY1]);
    tY(end) = mY1/norm([mX1 mY1]);
    Nend = [0 -1; 1 0]*[mX1;mY1];
    nX(end) = Nend(1)/norm(Nend);
    nY(end) = Nend(2)/norm(Nend);
    
    % columnwise
    X = X';
    Y = Y';
    tX = tX';
    tY = tY';
    nX = nX';
    nY = nY';
end