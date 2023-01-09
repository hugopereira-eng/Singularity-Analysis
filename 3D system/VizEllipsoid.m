%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-CMG cluster visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 54.73*pi/180;    % Pyramid skew angle

%% Ellipsoid - plots

% Roof array
% g = [pi/2;0;0;-pi/2];      % Rank-1
% g = [pi/2;0;0;0];          % Rank-2
% g = [0;0;0;0];             % Rank-3
% J = RoofJacobian(g,h0);

% Pyramid array
% g = [pi/2;pi;3*pi/2;0];      % Rank-2
g = [0;0;0;0];             % Rank-3
J = PyramidJacobian(g,h0,beta);

% SVD
[S,Sigma,V] = svd(J);

xc = 0;
yc = 0;
zc = 0;

xr = Sigma(1,1);
yr = Sigma(2,2);
zr = Sigma(3,3);
[X,Y,Z] = ellipsoid(xc,yc,zc,xr,yr,zr);
Xc = reshape(X,[],1);
Yc = reshape(Y,[],1);
Zc = reshape(Z,[],1);
coord = [Xc, Yc, Zc];
Rot = S*coord';

Xr = Rot(1,:);
Yr = Rot(2,:);
Zr = Rot(3,:);

Xs = reshape(Xr, 21,21);
Ys = reshape(Yr, 21,21);
Zs = reshape(Zr, 21,21);
C = sqrt(Xs.^2+Ys.^2+Zs.^2);

figure

mesh(Xs,Ys,Zs,C);

xlabel('$$x$$ axis','Interpreter','latex','FontSize',15)
ylabel('$$y$$ axis','Interpreter','latex','FontSize',15)
zlabel('$$z$$ axis','Interpreter','latex','FontSize',15)
title('Principal Components Visualization','FontSize',15)
% subtitle('Roof array','FontSize',12)
subtitle('Pyramid array','FontSize',12)
grid off
box off
axis square
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;

% view(0, 90)
