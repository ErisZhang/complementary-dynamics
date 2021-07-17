% V          vertex positions
% I          indices of each vertex's closest cluster
% new_C new  Catmull handles
% sample_num number of samples for Catmull interpolation
% Vr         deformed positions by wire deformer
% X Y        updated cluster positions

function [Vr,X,Y] = wire_deform(V,I,new_C,sample_num,X0,Y0,tX0,tY0,nX0,nY0)
    [X,Y,tX,tY,nX,nY] = catmull_rom_handle(new_C(:,1),new_C(:,2),sample_num);
    Vr = zeros(size(V));
    % update position
    for i = 1:size(V,1)
        v = V(i,:);
        idx = I(i);
        Vr(i,:) = [tX(idx) nX(idx); tY(idx) nY(idx)]*inv([tX0(idx) nX0(idx); tY0(idx) nY0(idx)])*(v-[X0(idx) Y0(idx)])'+[X(idx) Y(idx)]';
    end
end