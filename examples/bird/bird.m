clear;
warning('off','all');
addpath('../../matlab-include/utils');
set_project_path();

% read json
json = jsondecode(fileread('./bird_mouth_data/bird_mouth.json'));
cellfun(@(x,y) assignin('base',x,y),fieldnames(json),struct2cell(json));

%[V,F] = load_mesh(mesh_file);
[V,F,UV,TF] = readOBJ(mesh_file);
[C,~,P,BE,CE,~] = readTGF(handle_file);
V = V(:,1:2);
C = C(:,1:2);
[~,b] = farthest_points(V,5);
H = kharmonic(V,F,b,linspace(0,1,numel(b))');

% import handle transformations
T_list = read_2d_bone_anim(anim_file,C,BE);
[W] = save_bbw_bone('./bird_mouth_data/bird_mouth.dmat',V,F,C,P,BE);
VM = lbs_matrix(V,W);

clf;
hold on;
t = tsurf(F,V,'FaceColor',blue,'FaceAlpha',0.8,'EdgeAlpha',0.8);
hold off;
axis equal;
axis(5*[-0.5 2.3 -0.3 0.3]);
axis manual;
drawnow;

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

g = 0*vec(repmat([0 -9.8],size(V,1),1)); % didn't add any external gravity in this example
mask = M * phi;
Aeq = A'*mask;
Beq = zeros(size(A', 1),1);

% Bartels
[lambda, mu] = emu_to_lame(YM*ones(size(F,1),1), pr*ones(size(F,1),1));
dX = linear_tri2dmesh_dphi_dX(V,F);
areas = triangle_area(V,F);

% arap
energy_func = @(a,b,c,d,e) linear_tri2dmesh_arap_q(a,b,c,d,e,[0.5*lambda,mu]);
gradient_func = @(a,b,c,d,e) linear_tri2dmesh_arap_dq(a,b,c,d,e,[0.5*lambda,mu]);
hessian_func = @(a,b,c,d,e) linear_tri2dmesh_arap_dq2(a,b,c,d,e,[0.5*lambda,mu],'fixed');


% Don't move
U = vec(zeros(size(V)));
Ud = vec(zeros(size(V)));
Uc = vec(zeros(size(V)));

T = zeros(1,size(T_list,1));

for ai=1:size(T_list,1)
    
    tic;
    
    % read a new frame
    TCol = T_list{ai,1};
    UCol = VM*TCol;
    Ur = reshape(UCol, size(V));
    % normalize the mesh [0 1 0 1]: units are meters
    Ur = Ur-mesh_center;
    Ur = Ur/mesh_scale;
    
    Ur = vec(Ur-V);
     

    Ud0 = Ud;
    U0 = U;


      if with_dynamics


        % instead of one direct solve: do newton here
        max_iter = 20;
    %         Uc = vec(zeros(size(V))); % initial guess for Uc
        for i = 1 : max_iter
          alpha = 1;
          p = 0.5;
          c = 1e-8;

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

      U = Ur + Uc;
      Ud = (U-U0)/dt;
      
      t.Vertices = V+reshape(U,size(V,2),size(V,1))';
      title(sprintf('%d',ai),'Fontsize',20);
      drawnow;
end

