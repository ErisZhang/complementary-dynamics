function [T_list] = read_2d_bone_anim(filename,C,BE)
    
    fileID=fopen(filename,'r');
    nf=fscanf(fileID,'%d',1);
    m=fscanf(fileID,'%d',1);
    TM=fscanf(fileID,'%f', [3*m nf+1]);
    fclose(fileID);
    TM = TM';
    
    P = bone_parents(BE);
    
    d2r = pi/180;
    
    T_list = cell(nf,1);
    
    
    Tlr = zeros(m,2);
    for b=1:m
        Tlr(b,:) = TM(1,3*(b-1)+1:3*(b-1)+2);
    end

    function fk2d(b)
        if ~computed(b)
          p = P(b);
          if p < 1
            theta = TMi(:,3*(b-1)+3);
            %theta = TMi(b);
            Ti(3*(b-1)+1:3*b,:) = transform2d(theta*d2r,C(BE(b,1),1),C(BE(b,1),2));
          else
            fk2d(p);
            theta = TMi(:,3*(b-1)+3);
            %theta = TMi(b);
            Ti(3*(b-1)+1:3*b,:) = Ti(3*(p-1)+1:3*p,:)*transform2d(theta*d2r,C(BE(b,1),1),C(BE(b,1),2));
          end
          computed(b) = true;
        end
    end
    

    for ai=1:nf
        computed = false(m,1);
        TMi = TM(ai+1,:);
        Ti = zeros(3*m,3);
        for b =1:m
            fk2d(b);
        end
        T = zeros(2,3,m);
        for b =1:m
            if(P(b)<1)
                Tl = TM(ai+1,3*(b-1)+1:3*(b-1)+2)-Tlr(b,:);
                Ti(3*(b-1)+1:3*(b-1)+2,3) = Ti(3*(b-1)+1:3*(b-1)+2,3) + Tl';
            end
            T(:,:,b) = Ti(3*(b-1)+1:3*(b-1)+2,:);
        end
        T_list{ai,1} = reshape(permute(T,[3 1 2]),m*6,1);
    end

end

