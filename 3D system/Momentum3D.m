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
samples = 40;
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
%                 h(:,j) = RoofMomentum(g,h0);
                h(:,j) = PyramidMomentum(g,h0,beta);
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
% subtitle('Pyramid array','FontSize',12)
subtitle('$\beta = 54.73^{\circ}$','Interpreter','latex','FontSize',12);
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

%% Singular momentum space
n  = 10;
g1 = sym('g1', 'real');
g2 = sym('g2', 'real');
g3 = linspace(-pi,pi,n);
g4 = linspace(-pi,pi,n);
s = [];
for i3 = 1:length(g3)
    for i4 = 1:length(g4)   
        D1 = sin(beta)*(sin(g3(i3))*sin(g4(i4))*cos(g1) + cos(g3(i3))*sin(g4(i4))*sin(g1) + ...
            cos(beta)*(cos(g3(i3))*cos(g4(i4))*sin(g1) - sin(g3(i3))*cos(g4(i4))*cos(g1)) + ...
            2*cos(beta)^2*cos(g3(i3))*cos(g4(i4))*cos(g1));
        D2 = sin(beta)*(sin(g2)*sin(g3(i3))*cos(g4(i4)) + cos(g2)*sin(g3(i3))*sin(g4(i4)) + ...
            cos(beta)*(cos(g2)*cos(g3(i3))*sin(g4(i4)) - sin(g2)*cos(g3(i3))*cos(g4(i4))) + ...
            2*cos(beta)^2*cos(g2)*cos(g3(i3))*cos(g4(i4)));
        eq1 = D1 == 0;
        eq2 = D2 == 0;
        [solg1,solg2] = solve(eq1,eq2,g1,g2);
        s1 = [double(solg1(1)) double(solg1(2)) double(solg1(3)) double(solg1(4))];
        s2 = [double(solg2(1)) double(solg2(2)) double(solg1(3)) double(solg1(4))];
        s = [s [real(s1); real(s2)]];    
    end
end

%% Plot: singular momentum space
h_t = [];
g32 = [];
for i = 1:10
    for j = 1:10
        g32 = [g32 [g3(i); g4(j)]];
    end
end
for i = 1:100
    g = [s(1,4*i-3),s(2,4*i-3),g32(1,i),g32(2,i)];
    h = PyramidMomentum(g,h0,beta);
    h_t = [h_t h];
end
for i = 1:100
    g = [s(1,4*i-2),s(2,4*i-2),g32(1,i),g32(2,i)];
    h = PyramidMomentum(g,h0,beta);
    h_t = [h_t h];
end
C = real(sqrt(h_t(1,:).^2+h_t(2,:).^2+h_t(3,:).^2));
figure
scatter3(h_t(1,:),h_t(2,:),h_t(3,:), 10, C, 'fill')
c = colorbar;
colormap default;
c.Label.String = 'Magnitude |h|';
xlabel('$\dot{h}_x$ [Nms]','Interpreter','latex','FontSize',15);
ylabel('$\dot{h}_y$ [Nms]','Interpreter','latex','FontSize',15);
zlabel('$\dot{h}_z$ [Nms]','Interpreter','latex','FontSize',15);
title('Singular momentum space','Interpreter','latex','FontSize',15);
subtitle('Pyramid array','FontSize',12);
xlim([-pi pi])
ylim([-pi pi])
zlim([-pi pi])
grid off
axis equal