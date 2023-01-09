%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-CMG momentum visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 54.73*pi/180;    % Pyramid skew angle

%% Angular momentum - 3D

samples = 10;
g1 = linspace(-pi,pi,samples);
g2 = linspace(-pi,pi,samples);
g3 = linspace(-pi,pi,samples);
g4 = linspace(-pi,pi,samples);     

j = 1;
h = zeros(3,samples^4);
for i1 = 1:length(g1)
    for i2 = 1:length(g2)
        for i3 = 1:length(g3)
            for i4 = 1:length(g4)
                g = [g1(i1);g2(i2);g3(i3);g4(i4)];
                h(:,j) = RoofMomentum(g,h0);
%                 h(:,j) = PyramidMomentum(g,h0,beta);
                j = j + 1;
            end
        end
    end
end
L = sqrt(length(h));
xm = reshape(h(1,:),L,L);
ym = reshape(h(2,:),L,L);
zm = reshape(h(3,:),L,L);
C = sqrt(xm.^2+ym.^2+zm.^2);
% mesh(xm,ym,zm,C)
scatter3(h(1,:),h(2,:),h(3,:),10,reshape(C,1,L^2),'fill')
xlabel('$$h_x/h_0$$','Interpreter','latex','FontSize',15)
ylabel('$$h_y/h_0$$','Interpreter','latex','FontSize',15)
zlabel('$$h_z/h_0$$','Interpreter','latex','FontSize',15)
title('Angular momentum envelope','FontSize',15)
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