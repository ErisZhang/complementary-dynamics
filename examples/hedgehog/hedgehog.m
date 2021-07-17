clear all;
warning('off','all');
addpath('../../matlab-include/utils');
set_project_path();

% read json
json = jsondecode(fileread('./hedgehog_data/hedgehog_data.json'));
cellfun(@(x,y) assignin('base',x,y),fieldnames(json),struct2cell(json));

has_user_input = true;

[V,F] = load_mesh(mesh_file);
V = V(:,1:2);
BE = [];
CE = [];
PE = [];
I = readDMAT(constraint_file);
C = V(I,:);
P = 1:size(C,1); % default value for point handle


% normalize the mesh
Vcenter = (min(V)+max(V))*0.5;
V = V-Vcenter;
Vbound = max(V(:));
V = V/Vbound;
C = C-Vcenter;
C = C/Vbound;
[~,b] = farthest_points(V,5);
H = kharmonic(V,F,b,linspace(0,1,numel(b))');

% % Compute boundary conditions
[b,bc] = boundary_conditions(V,F,C,P,BE,CE);
% Compute weights
W = biharmonic_bounded(V,F,b,bc,'OptType','quad');
% Normalize weights
W = W./repmat(sum(W,2),1,size(W,2));


vec = @(X) reshape(X',size(X,1)*size(X,2),1);
M = massmatrix_xyz(V,F);

A = lbs_matrix_xyz(V,W);

% no external force use the default D matrix
[phi,Em] = default_D_matrix(V,F);


clf;
hold on;
t = tsurf(F,V,'FaceColor',blue,'FaceAlpha',0.8,'EdgeAlpha',0.8);
% visualize the handles
C_plot = scatter3( ...
C(:,1),C(:,2),0*C(:,1), ... 
'o','MarkerFaceColor',[1 1 1], 'MarkerEdgeColor','k',...
'LineWidth',2,'SizeData',50);
axis equal;
expand_axis(3);
% axis(5*[-0.5 1.5 -0.3 0.3]);
axis manual;
drawnow;

g = 0*vec(repmat([0 -9.8],size(V,1),1)); % no gravity

mask = M * phi;
Aeq = A'*mask;
Beq = zeros(size(A', 1),1);

% Bartels
[lambda, mu] = emu_to_lame(YM*ones(size(F,1),1), pr*ones(size(F,1),1));
dX = linear_tri2dmesh_dphi_dX(V,F);
areas = triangle_area(V,F);

energy_func = @(a,b,c,d,e) linear_tri2dmesh_arap_q(a,b,c,d,e,[0.5*lambda,mu]);
gradient_func = @(a,b,c,d,e) linear_tri2dmesh_arap_dq(a,b,c,d,e,[0.5*lambda,mu]);
hessian_func = @(a,b,c,d,e) linear_tri2dmesh_arap_dq2(a,b,c,d,e,[0.5*lambda,mu],'fixed');


% Don't move
U = vec(zeros(size(V)));
Ud = vec(zeros(size(V)));
Uc = vec(zeros(size(V)));

g = 0 * vec(repmat([9.8 0],size(V,1),1));

user_T = readDMAT(anim_file);

with_dynamics = true;

frame_num = size(user_T, 1);
T = zeros(1,frame_num);

% while loop for animation
for iter = 1:frame_num
    
  Ud0 = Ud;
  U0 = U;
  
  Z = user_T(iter,:)';
  new_Vr = reshape(A*Z,size(V,2),size(V,1))';
  new_C = new_Vr(I,:);
  
  Ur = A*Z-vec(V); % rig displacement

  if with_dynamics

    % instead of one direct solve: do newton here
    max_iter = 20;
%     Uc = vec(zeros(size(V))); % initial guess for Uc
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
      if tmp_g'*dUc > -1e-6
         break;
      end

      % backtracking line search
      alpha = newton_line_search(f,tmp_g,dUc,Ur,Uc);
      Uc = Uc + alpha * dUc;
   
    end
    
  end
  
  U = Ur + with_dynamics * Uc;
  Ud = (U-U0)/dt;
  
  t.Vertices = V+reshape(U,size(V,2),size(V,1))';
  
  set(C_plot,'XData',new_C(:,1));
  set(C_plot,'YData',new_C(:,2));

  title(sprintf('%d',iter),'Fontsize',20);
  drawnow;
  
%     %%%%%%%%%%%%%%% save file %%%%%%%%%%%%%%%
%     U_output = V + reshape(U,size(V,2),size(V,1))';
%     file_name = "./output/";
% 
%     saveOBJ(U_output,F,file_name,iter);
%     %%%%%%%%%%%%%%% save file %%%%%%%%%%%%%%%

end


