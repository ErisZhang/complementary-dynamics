clear all;
warning('off','all');
addpath('../../matlab-include/utils');
set_project_path();

% read json
json = jsondecode(fileread('./walrus_data/walrus.json'));
cellfun(@(x,y) assignin('base',x,y),fieldnames(json),struct2cell(json));

[V,F] = load_mesh(mesh_file);
V = V(:,1:2);
[C,E] = readTGF(skeleton_file);
C = C(:,1:2);


% normalize vertices positions to be [-1 1 -1 1]
Vcenter = (min(V)+max(V))*0.5;
V = V-Vcenter; Vbound = max(V(:)); V = V/Vbound;
C = C-Vcenter; C = C/Vbound;

CE = E; % CE and E are the same in this case
P = 1:size(C,1);
BE = [];

% compute Harmonic Coordinates W
[TV,TF] = triangle([V;C],size(V,1) + CE,[],'Flags','-q31a0.01');
[b,bc] = boundary_conditions(TV,TF,C,P,[],CE);
TW = kharmonic(TV,TF,b,bc,1);
W = TW(1:size(V,1),:);

vec = @(X) reshape(X',size(X,1)*size(X,2),1);
M = massmatrix_xyz(V,F);

% harmonic coordinates
A = zeros(size(W,1)*2,size(W,2)*2);
for i = 1:size(W,1)
    A(2*(i-1)+1:2*(i-1)+2,:) = repdiag(W(i,:),2);
end

% no external force use the default D matrix
[phi,Em] = default_D_matrix(V,F);

% viewr
clf;
hold on;
t = tsurf(F,V,'FaceColor',blue,'FaceAlpha',0.8,'EdgeAlpha',0.8);
C_plot = scatter3( ...
C(:,1),C(:,2),0*C(:,1), ... 
'o','MarkerFaceColor',[255 165 0]/255, 'MarkerEdgeColor',[255 165 0]/255,...
'LineWidth',1,'SizeData',50);
% draw the cage
L_plot = line('XData',C([CE(:,1);1],1),'YData',C([CE(:,1);1],2),'LineWidth',2,'Color',[255 165 0]/255);
hold off;
axis equal;
expand_axis(1.5);
axis manual;
drawnow;


mask = M * phi;
Aeq = A' * mask;
Beq = zeros(size(A', 1),1);

% Bartels
[lambda, mu] = emu_to_lame(YM*ones(size(F,1),1), pr*ones(size(F,1),1));
dX = linear_tri2dmesh_dphi_dX(V,F);
areas = triangle_area(V,F);
energy_func = @(a,b, c, d, e) linear_tri2dmesh_neohookean_q(a,b,c,d,e,[0.5*mu, 0.5*lambda]);
gradient_func = @(a,b, c, d, e) linear_tri2dmesh_neohookean_dq(a,b,c,d,e,[0.5*mu, 0.5*lambda]);
hessian_func = @(a,b, c, d, e) linear_tri2dmesh_neohookean_dq2(a,b,c,d,e,[0.5*mu, 0.5*lambda], 'fixed');


% Don't move
U = vec(zeros(size(V)));
Ud = vec(zeros(size(V)));
Uc = vec(zeros(size(V)));

g = 0 * vec(repmat([9.8 0],size(V,1),1)); % no gravity

user_C = readDMAT(anim_file);

frame_num = size(user_C,1);
T = zeros(1,frame_num);

% while loop for animation
for iter = 1:frame_num
    
  Ud0 = Ud;
  U0 = U;
  
  new_C = user_C(iter,:)';
  
  Ur = A * new_C  - vec(V); % rig displacement

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


      tmp_H = 0.5 * (tmp_H + tmp_H');
      dUc = speye(size(tmp_H,1),size(tmp_H,1)+size(Aeq,1)) * ([tmp_H Aeq';Aeq sparse(size(Aeq,1),size(Aeq,1))] \ [-tmp_g;Beq]);

      if tmp_g' * dUc > -1e-6
          break;
      end
      
      % backtracking line search
      alpha = newton_line_search(f,tmp_g,dUc,Ur,Uc);
      Uc = Uc + alpha * dUc;
      
    end
  end
  
  U = Ur + with_dynamics * Uc;
  Ud = (U-U0)/dt;
  
  % update the visualization
  set(t,'Vertices',V+reshape(U,size(V,2),size(V,1))');
  new_C = reshape(new_C,size(new_C,1)/2,2);
  set(C_plot,'XData',new_C(:,1));
  set(C_plot,'YData',new_C(:,2));
  set(L_plot,'XData',new_C([CE(:,1);1],1));
  set(L_plot,'YData',new_C([CE(:,1);1],2));
  
  title(sprintf('%d',iter),'Fontsize',20);
  drawnow;

end

