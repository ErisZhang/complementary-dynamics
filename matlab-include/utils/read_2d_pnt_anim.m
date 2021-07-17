function [T_list] = read_2d_pnt_anim(filename,C,PI)


    %TM = readMatrix(filename);
    
    fileID=fopen(filename,'r');
    nf=fscanf(fileID,'%d',1);
    m=fscanf(fileID,'%d',1);
    TM=fscanf(fileID,'%f', [3*m nf+1]);
    fclose(fileID);
    TM = TM';

    d2r = pi/180;
    
    T_list = cell(nf,1);


    Tlr = zeros(m,2);
    for b=1:m
        Tlr(b,:) = TM(1,3*(b-1)+1:3*(b-1)+2);
    end
    

    for ai=1:nf
        T = zeros(2,3,m);
        for b =1:m
            Tl = TM(ai+1,3*(b-1)+1:3*(b-1)+2)-Tlr(b,:);
            theta = TM(ai+1,3*(b-1)+3);
            TR = transform2d(theta*d2r,C(PI(b),1),C(PI(b),2));
            TR(1:2,3) = TR(1:2,3) + Tl';
            T(:,:,b) = TR(1:2,:);
        end
        T_list{ai,1} = reshape(permute(T,[3 1 2]),m*6,1);
    end

end

