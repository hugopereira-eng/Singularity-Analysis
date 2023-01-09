%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-CMG torque visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 54.73*pi/180;    % Pyramid skew angle

%% Torque envelope - 2D

% Roof array
% g = [pi/2;0;0;-pi/2];      % Rank-1
% g = [pi/2;0;0;0];          % Rank-2
% g = [0;0;0;0];             % Rank-3
% J = RoofJacobian(g,h0);

% Pyramid array
% g = [pi/2;pi;3*pi/2;0];      % Rank-2
g = [0;0;0;0];             % Rank-3
J = PyramidJacobian(g,h0,beta);

samples = 10;
gr1 = linspace(-1.5,1.5,samples);
gr2 = linspace(-1.5,1.5,samples);
gr3 = linspace(-1.5,1.5,samples);
gr4 = linspace(-1.5,1.5,samples);

hDot = zeros(3,samples^4);
j = 1;
for i1 = 1:length(gr1)
    for i2 = 1:length(gr2)
        for i3 = 1:length(gr3)
            for i4 = 1:length(gr4)
                gr = [gr1(i1);gr2(i2);gr3(i3);gr4(i4)];
                hDot(:,j) = J*gr;
                j = j + 1;
            end
        end
    end
end

L = sqrt(length(hDot));
C = sqrt(hDot(1,:).^2+hDot(2,:).^2+hDot(3,:).^2);
scatter3(hDot(1,:),hDot(2,:),hDot(3,:),10,reshape(C,1,L^2),'fill')
hxl = xlabel('$$\dot{h}_x/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15);
hxl.Position=[4 -1 -8];
hyl = ylabel('$$\dot{h}_y/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15);
hyl.Position=[-5 2.5 -5.5];
zlabel('$$\dot{h}_z/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15)
title('Torque envelope','FontSize',15)
% subtitle('Roof array','FontSize',12)
subtitle('Pyramid array','FontSize',12)
grid off
box off
axis square
xlim([-4 4])
ylim([-4 4])
zlim([-4 4])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;
% view(0,0)
% view(0,90)