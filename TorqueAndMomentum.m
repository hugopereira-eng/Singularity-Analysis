%% Torque envelope - 2D

h0 = 1;
g = [0;2*pi/3;-2*pi/3];
% g = [0;0;0];
alpha = 30*pi/180;

J = compute_Jac_triangle(g,h0,alpha);

points = 100;
gr1 = linspace(-1.5,1.5,points);
gr2 = linspace(-1.5,1.5,points);
gr3 = linspace(-1.5,1.5,points);

j = 1;
for i1 = 1:length(gr1)
    for i2 = 1:length(gr2)
        for i3 = 1:length(gr3)
            gr = [gr1(i1);gr2(i2);gr3(i3)];
            h_dot(:,j) = J*gr;
            j = j + 1;
        end
    end
end

C = sqrt(h_dot(1,:).^2+h_dot(2,:).^2);
scatter(h_dot(1,:),h_dot(2,:),10,C,'filled')
xlabel('$$\dot{h}_x/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15)
ylabel('$$\dot{h}_y/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15)
title('Torque envelope','FontSize',15)
subtitle('Triangular array','FontSize',12)
grid off
box off
axis square
xlim([-3.5 3.5])
ylim([-3.5 3.5])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;

%% Torque envelope - 3D

h0 = 1;
% g = [pi/2;0;0;-pi/2];
g = [pi/2;0;0;0];
% g = [0;0;0;0];

J = compute_Jac_roof(g,h0);

points = 60;
gr1 = linspace(-1.5,1.5,points);
gr2 = linspace(-1.5,1.5,points);
gr3 = linspace(-1.5,1.5,points);
gr4 = linspace(-1.5,1.5,points);

j = 1;
for i1 = 1:length(gr1)
    for i2 = 1:length(gr2)
        for i3 = 1:length(gr3)
            for i4 = 1:length(gr4)
                gr = [gr1(i1);gr2(i2);gr3(i3);gr4(i4)];
                h_dot(:,j) = J*gr;
                j = j + 1;
            end
        end
    end
end

L = sqrt(length(h_dot));
% xm = reshape(h_dot(1,:),L,L);
% ym = reshape(h_dot(2,:),L,L);
% zm = reshape(h_dot(3,:),L,L);
% C = sqrt(xm.^2+ym.^2+zm.^2);
% mesh(xm,ym,zm,C)
C = sqrt(h_dot(1,:).^2+h_dot(2,:).^2+h_dot(3,:).^2);
scatter3(h_dot(1,:),h_dot(2,:),h_dot(3,:),10,reshape(C,1,L^2),'fill')
hxl = xlabel('$$\dot{h}_x/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15);
hxl.Position=[4 -1 -8];
hyl = ylabel('$$\dot{h}_y/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15);
hyl.Position=[-5 2.5 -5.5];
zlabel('$$\dot{h}_z/h_0$$ [$$s^{-1}$$]','Interpreter','latex','FontSize',15)
title('Torque envelope','FontSize',15)
subtitle('Roof array','FontSize',12)
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

%% Angular momentum - 2D

points = 100;
g1 = linspace(-pi,pi,points);
g2 = linspace(-pi,pi,points);
g3 = linspace(-pi,pi,points);

alpha = 30*pi/180;
h0 = 1;     

j = 1;
h = zeros(2,points^3);
for i1 = 1:length(g1)
    for i2 = 1:length(g2)
        for i3 = 1:length(g3)
            g = [g1(i1);g2(i2);g3(i3)];
            h(:,j) = compute_h_triangle(h0,alpha,g);
            j = j + 1;
        end
    end
end
C = sqrt(h(1,:).^2+h(2,:).^2);
scatter(h(1,:),h(2,:),[],C,'filled')
xlabel('$$h_x/h_0$$','Interpreter','latex','FontSize',15)
ylabel('$$h_y/h_0$$','Interpreter','latex','FontSize',15)
title('Angular momentum envelope','FontSize',15)
subtitle('Triangular array','FontSize',12)
grid off
box off
axis square
xlim([-3.5 3.5])
ylim([-3.5 3.5])
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 12;

%% Angular momentum - 3D

points = 100;
g1 = linspace(-pi,pi,points);
g2 = linspace(-pi,pi,points);
g3 = linspace(-pi,pi,points);
g4 = linspace(-pi,pi,points);

h0 = 1;      

j = 1;
h = zeros(3,points^4);
for i1 = 1:length(g1)
    for i2 = 1:length(g2)
        for i3 = 1:length(g3)
            for i4 = 1:length(g4)
                g = [g1(i1);g2(i2);g3(i3);g4(i4)];
                h(:,j) = compute_h_roof(h0,g);
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
mesh(xm,ym,zm,C)
% scatter3(h(1,:),h(2,:),h(3,:),10,reshape(C,1,L^2),'fill')
xlabel('$$h_x/h_0$$','Interpreter','latex','FontSize',15)
ylabel('$$h_y/h_0$$','Interpreter','latex','FontSize',15)
zlabel('$$h_z/h_0$$','Interpreter','latex','FontSize',15)
title('Angular momentum envelope','FontSize',15)
subtitle('Roof array','FontSize',12)
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
%% Functions

function momentum = compute_h_triangle(h0,alpha,g)
    momentum = h0*[cos(alpha+g(1))+sin(g(2))-cos(alpha-g(3));
                   sin(alpha+g(1))-cos(g(2))+sin(alpha-g(3))];
end

function momentum = compute_h_roof(h0,g)
    momentum = h0*[sin(g(1))+cos(g(2));
                   -cos(g(3))+sin(g(4));
                   -cos(g(1))+sin(g(2))+sin(g(3))+cos(g(4))];
end

function J = compute_Jac_triangle(g,h0,alpha)
J = h0*[-sin(alpha + g(1))   cos(g(2))  -sin(alpha - g(3));
         cos(alpha + g(1))   sin(g(2))  -cos(alpha - g(3))];
end

function J = compute_Jac_roof(g,h0)
J = h0*[ cos(g(1))   -sin(g(2))             0            0;
                 0            0     sin(g(3))    cos(g(4));
         sin(g(1))    cos(g(2))     cos(g(3))   -sin(g(4))];
end
