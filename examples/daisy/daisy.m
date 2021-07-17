clear;
warning('off','all');
addpath('../../matlab-include/utils');
set_project_path();

% read json
json = jsondecode(fileread('./daisy_data/daisy.json'));
cellfun(@(x,y) assignin('base',x,y),fieldnames(json),struct2cell(json));


[V,F] = load_mesh(mesh_file);
[C,~,PI,BE,CE,~] = readTGF(handle_file);
V(:,3)=[];
C(:,3)=[];
[~,b] = farthest_points(V,5);
H = kharmonic(V,F,b,linspace(0,1,numel(b))');

T_list = read_2d_pnt_anim(anim_file,C,PI);
W = save_bbw_pnt(weight_file, V,F,C);

I = eye(size(PI,2));
CM = lbs_matrix(C,I);
VM = lbs_matrix(V,W);


% normalize the mesh
mesh_center = (max(V)+min(V))/2;
V = V - mesh_center;
mesh_scale = max(V(:));
V = V/mesh_scale;
C = C - mesh_center;
C = C/mesh_scale;

vec = @(X) reshape(X',size(X,1)*size(X,2),1);

M = massmatrix_xyz(V,F);
A = lbs_matrix_xyz(V,W);

% no external force use the default D matrix
[phi,Em] = default_D_matrix(V,F);

mask = M * phi;

g = 0*vec(repmat([0 -9.8],size(V,1),1)); % didn't add any external gravity in this example
Beq = zeros(size(A', 1),1);

% Bartels
[lambda, mu] = emu_to_lame(YM*ones(size(F,1),1), pr*ones(size(F,1),1));
dX = linear_tri2dmesh_dphi_dX(V,F);
areas = triangle_area(V,F);

energy_func = @(a,b,c,d,e) linear_tri2dmesh_arap_q(a,b,c,d,e,[0.5*lambda,mu]);
gradient_func = @(a,b,c,d,e) linear_tri2dmesh_arap_dq(a,b,c,d,e,[0.5*lambda,mu]);
hessian_func = @(a,b,c,d,e) linear_tri2dmesh_arap_dq2(a,b,c,d,e,[0.5*lambda,mu],'fixed');

% plot the wire
sample_num = 30;
[X0,Y0,tX0,tY0,nX0,nY0] = catmull_rom_handle(C(:,1),C(:,2),sample_num); % reference frame
Wv = [X0 Y0];
[I] = cluster_pnts_to_closest(V,Wv);

% Don't move
U = vec(zeros(size(V)));
Ud = vec(zeros(size(V)));
Uc = vec(zeros(size(V)));


clf;
hold on;
t = tsurf(F,V,'FaceColor',blue,'FaceAlpha',0.8,'EdgeAlpha',0.8);
C_plot = scatter3( ...
C(:,1),C(:,2),0*C(:,1), ... 
'o','MarkerFaceColor',[1 1 1], 'MarkerEdgeColor','k',...
'LineWidth',2,'SizeData',50);
L_plot = line('XData',X0,'YData',Y0,'LineWidth',1.5,'Color',[1 1 1]);
hold off;
axis equal;
expand_axis(3);
% axis(5*[-0.5 1.5 -0.3 0.3]);
axis manual;
drawnow;


avg_edge_len = avgedge(V,F); % before loop

Ti = zeros(1,size(T_list,1));

for ai = 1:size(T_list,1)
    
    % read a new frame
    TCol = T_list{ai,1};
    
    CCol = CM*TCol;
    new_C = reshape(CCol, size(C));
    new_C = new_C-mesh_center;
    new_C = new_C/mesh_scale;
    
    [Vr,X,Y] = wire_deform(V,I,new_C,sample_num,X0,Y0,tX0,tY0,nX0,nY0);
    Ud0 = Ud;
    U0 = U;
    Uc0 = Uc;
    Ur = vec(Vr)-vec(V); % rig displacement


      if with_dynamics
      
        % compute the jacobian
        h = 1e-8;
        J = zeros(size(V,1)*size(V,2), size(new_C,1)*size(new_C,2));
        Vtmp = zeros(size(V));
        % finite difference
        for i = 1:size(J,2)
            row = ceil(i/2);
            col = (i-row*2)+2;
            new_C(row,col) = new_C(row,col) + h;
            [Vtmp,~,~] = wire_deform(V,I,new_C,sample_num,X0,Y0,tX0,tY0,nX0,nY0);
            J(:,i) = (vec(Vtmp) - vec(Vr))/h;
            new_C(row,col) = new_C(row,col) - h;
        end
        Aeq = J' * mask;
        Beq = zeros(size(J',1),1);

        
        % instead of one direct solve: do newton here
        max_iter = 20;
%         Uc = vec(zeros(size(V))); % initial guess for Uc
        for i = 1 : max_iter

          % total energy = gravitational potential energy + kinetic energy +
          % elastic potential energy
          G = gradient_func(V,F,vec(V)+Ur+Uc,dX,areas);
          K = hessian_func(V,F,vec(V)+Ur+Uc,dX,areas);
    
          f = @(Ur,Uc) -(M*g)'*(vec(V)+Ur+Uc) + 0.5*(Ur+Uc-U0-dt*Ud0)'*M/(dt*dt)*(Ur+Uc-U0-dt*Ud0) + ...
              energy_func(V,F,vec(V)+Ur+Uc,dX,areas);
            
            
          tmp_H = M/(dt^2) + K;
          tmp_g = M/(dt^2) * (Ur+Uc) - M*(g+U0/dt^2+Ud0/dt) + G;

          % solve for the update
          tmp_H = 0.5 * (tmp_H + tmp_H');
          dUc = speye(size(tmp_H,1),size(tmp_H,1)+size(Aeq,1)) * ([tmp_H Aeq';Aeq sparse(size(Aeq,1),size(Aeq,1))] \ [-tmp_g;Beq]);

          % check for newton convergence criterian
          if tmp_g' * dUc > -1e-6
              break;
          end

          % backtracking line search
          alpha = newton_line_search(f,tmp_g,dUc,Ur,Uc);
          Uc = Uc + alpha * dUc;

        end
      end
      
      U = Ur + Uc;
      Ud = (U-U0)/dt;

      t.Vertices = V+reshape(U,size(V,2),size(V,1))';
      set(C_plot,'XData',new_C(:,1));
      set(C_plot,'YData',new_C(:,2));
      set(L_plot,'XData',X);
      set(L_plot,'YData',Y);
      title(sprintf('%d',ai),'Fontsize',20);
      drawnow;
     

end


