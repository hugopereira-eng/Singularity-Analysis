%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-CMG cluster visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
h0 = 1;                  % Angular momentum of each CMG
lambda = 30*pi/180;      % Triangle inner angle

%% Ellipse - plots
g = [0;2*pi/3;-2*pi/3];     % Non-singular configuration
% g = [0;0;0];                % Singularity
J = TriangleJacobianExtended(g,h0,lambda);

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

mesh(Xs,Ys,Zs,C);

xlabel('$$x$$-axis','Interpreter','latex','FontSize',15)
ylabel('$$y$$-axis','Interpreter','latex','FontSize',15)
zlabel('$$z$$-axis','Interpreter','latex','FontSize',15)
title('Principal Components Visualization','FontSize',15)
subtitle('Triangular array','FontSize',12)
grid off
box off
axis square
xlim([-1.5 1.5])
ylim([-1.5 1.5])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;

view(0, 90)

%% Function: Extended triangle Jacobian
function J = TriangleJacobianExtended(g,h0,lambda)
J = h0*[-sin(lambda + g(1))   cos(g(2))  -sin(lambda - g(3))  0;
         cos(lambda + g(1))   sin(g(2))  -cos(lambda - g(3))  0;
                        0           0                  0  0];
end
