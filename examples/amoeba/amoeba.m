clear all;
warning('off','all');
addpath('../../matlab-include/utils');
set_project_path();


% read json
json = jsondecode(fileread('./amoeba_data/amoeba.json'));
cellfun(@(x,y) assignin('base',x,y),fieldnames(json),struct2cell(json));


[TT,~,TA] = imread(image_file);
TT = im2double(TT);
[V,F] = bwmesh(TA,'SmoothingIters',10,'Tol',1);
UV = V;


% rendering setup
[~,b] = farthest_points(V,5);
H = kharmonic(V,F,b,linspace(0,1,numel(b))');
if rendering
    % mode = 'view-in-figure-2';
    % mode = 'gif';
    mode = 'png-sequence';
    switch mode
        case 'view-in-figure-2'
          close all;
          figure('WindowSTyle','docked');
          figure('WindowSTyle','docked');
          figure(1);
    end
end

% normalize the mesh
V = V(:,1:2);
V = V-(min(V)+max(V))*0.5;
V = V/max(V(:));

% linear elasticity
[~,~,data] = linear_elasticity(V,F,[],[],'Mu',mu,'Lambda',lambda);
K = data.K;
M = data.M;
vec = @(X) X(:);
g = 0*vec(repmat([0 -9.8],size(V,1),1)); % external force
VH = [V ones(size(V,1),1)];
A = sparse(repdiag(VH,size(V,2))); % single handle


% viewer
clf; 
hold on;
[X,Y] = meshgrid((min(V(:,1))-1):0.1:(max(V(:,1))+1),(min(V(:,2))-1):0.1:(max(V(:,2))+1));
[GF,GV] = surf2patch(X,Y,0*X);
q = quiver(GV(:,1),GV(:,2),0*GV(:,1),0*GV(:,2),'linewidth',2,'color',[145, 224, 255]/255);
t = tsurf(F,V,'EdgeColor','none');
hold off;
axis equal;
expand_axis(0.9);
axis manual;
drawnow;
% Remove the axis and make the background color white
set(gca,'Position',[0 0 1 1],'Visible','off');
set(gcf,'Color','w');


Z = zeros(size(A,2),1);
U = vec(zeros(size(V)));
Ud = vec(zeros(size(V)));
Uc = vec(zeros(size(V)));
T = [0 0];


for iter = 1:frame_num

  Z0 = Z;
  Ud0 = Ud;
  U0 = U;
  
  % script the input motion
  time = 0.025 * iter;
  T(end) = (sin(0.5*time))*0.6; 
  th = sin(0.5*time)*0.05;
  Z = vec([[cos(th) sin(th);-sin(th) cos(th)] - eye(2,2);T]);
  gfun = @(V) [max(sin(-6*time+V(:,1)*5),0) zeros(size(V,1),1)];
  g = vec(6e1*gfun(V));
  
  mqwf = {};
  mqwf.force_Aeq_li = true;
  [Uc,mqwf] = min_quad_with_fixed( ...
    0.5*(K+M/(dt^2)), ...
    M*(A*(Z - Z0))/(dt^2)-M*(g+Uc/dt^2+Ud0/dt), ...
    [],[], ...
    A' * M,zeros(size(A,2),1), ...
    mqwf);

  U = A*Z + Uc;
  Ud = (U-U0)/dt;
  
  % visualization
  t.Vertices = [V+reshape(U,size(V)) H];
  gGV = gfun(GV);
  set(q,'XData',GV(:,1),'YData',GV(:,2),'UData',gGV(:,1),'VData',gGV(:,2));
  title(sprintf('%d',iter),'Fontsize',20);
  drawnow;

  if rendering
      [IO,AA] = apply_texture_map(t,UV,TT);
      % show image or write to file
      switch mode
      case 'view-in-figure-2'
        figure(2);
        imshow(IO);
        set(gca,'Position',[0 0 1 1]);
        drawnow;
        figure(1);
      case 'gif'
        filename = './output/amoeba.gif';
        f = exist(filename,'file');
        [SIf,cm] = rgb2ind(IO,256);
        if ~f
          imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
        else
          imwrite(SIf,cm,filename,'WriteMode','append','Delay',0);
        end
      case 'png-sequence'
        imwrite(IO,sprintf('./output/amoeba%04d.png',iter),'Alpha',1*AA);
      end
  end
end
 

