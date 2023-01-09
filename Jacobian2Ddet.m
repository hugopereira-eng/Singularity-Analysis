close all
clear
clc

%% Parameters

h0 = 1;
beta = 30*pi/180;

%% Symbolic determinant expression

g1 = sym('g1', 'real');
g2 = sym('g2', 'real');
g3 = sym('g3', 'real');
g = [g1,g2,g3];  
J = computeJac(g,h0,beta );
% disp(simplify(det(J*J')))
r = [J(:,2) J(:,3)];
disp(simplify(det(r)^2))

%% Try values

g = [60; 0; 120]*pi/180;
J = computeJac(g,h0,beta);
fprintf("Determinant: %d \n",det(J*J'));

%% SVD

g = [0; 2*pi/3; -2*pi/3];
J = computeJac(g,h0,beta);
disp(J);
[U,S,V] = svd(J);

%% Gamma

gp = [];
g3 = -pi:0.001:pi;
for i = 1:length(g3)
    for j = 1:4
        if (j == 1)
            g1 = wrapToPi(g3(i) - pi/3);
            g2 = wrapToPi(g3(i) + pi/3);
        elseif (j == 2)
            g1 = wrapToPi(g3(i) - pi/3);
            g2 = wrapToPi(g3(i) + pi/3 + pi);
        elseif (j == 3)
            g1 = wrapToPi(g3(i) - pi/3 + pi);
            g2 = wrapToPi(g3(i) + pi/3);
        else
            g1 = wrapToPi(g3(i) - pi/3 + pi);
            g2 = wrapToPi(g3(i) + pi/3 + pi);
        end
        gp = [gp [g1;g2;g3(i)]];
    end
end

%% Angular momentum

h1h = []; index1h = [];
h3h = []; index3h = [];
for i = 1:length(gp)
    aux = MomentumCompute(gp(:,i),h0,beta);
    if round(norm(aux)) == h0
        h1h = [h1h aux];
        index1h = [index1h i];
    else 
        h3h = [h3h aux];
        index3h = [index3h i];
    end
end


%% Plot - Gamma plot

figure
scatter3(gp(1,index1h),gp(2,index1h),gp(3,index1h), 10,'fill','b')
hold on
scatter3(gp(1,index3h),gp(2,index3h),gp(3,index3h), 10,'fill','g')
hold on
plotcube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
title('Singularities envelope','Interpreter','latex','FontSize',15);
subtitle('Triangular array','FontSize',12);
legend('1h singularities','3h singularities','interpreter','latex','Location','southoutside','NumColumns',2);
axis equal
grid off
box on
view(335,20)

%% Plane x+y+z = 0

hold on
planeIndex = [];
for i = 1:length(gp)
    if (round(gp(1,i) + gp(2,i) + gp(3,i),2) == 0)
        planeIndex = [planeIndex gp(:,i)];
    end
end
scatter3(planeIndex(1,:),planeIndex(2,:),planeIndex(3,:), 15,'fill')

%% Plot - Singular angular momentum

figure
scatter(h1h(1,:),h1h(2,:), 10,'fill','b')
hold on
scatter(h3h(1,:),h3h(2,:), 10,'fill','g')
xlabel('$$h_x/h_0$$','Interpreter','latex','FontSize',15)
ylabel('$$h_y/h_0$$','Interpreter','latex','FontSize',15)
title('Singular momentum space','Interpreter','latex','FontSize',15);
subtitle('Triangular array','FontSize',12);
legend('1h singularities','3h singularities','interpreter','latex','Location','NorthWest');
axis square
xlim([-3.5 3.5])
ylim([-3.5 3.5])
grid off
box off

%% Optimal energy solution

% torque = [6;6];
% v1 = [-2*pi/3,-pi,-pi/3];
% v2 = [2*pi/3,pi/3,pi];
% pt = [2*pi/3,2*pi/3,pi];
% distance = point_to_line(pt,v1,v2);
% J = computeJac(pt,h0,beta);
% gDot = J'/(J*J')*torque;
% disp(gDot);
% disp(norm(gDot));
% disp(distance);

%% Singularity

% Desired torque
% torque = [sqrt(18);sqrt(18)];        
torque = [3;3];  
% 1st singular space
v11 = [-pi,-pi/3,-2*pi/3];
v12 = [pi/3,pi,2*pi/3];
% 2nd singular space
v21 = [-2*pi/3,-pi,-pi/3];
v22 = [2*pi/3,pi/3,pi];
% Straight discretization
indexvec = linspace(-0.5,1.5,1000);
% Allocation
d1       = zeros(1,length(indexvec));
d2       = zeros(1,length(indexvec));
gimbal   = zeros(3,length(indexvec));
gDot     = zeros(3,length(indexvec));
Norm2    = zeros(1,length(indexvec));
NormInf  = zeros(1,length(indexvec));

for i = 1:length(indexvec)
    pt = [-pi/3;pi/3;0] + indexvec(i)*[pi/3;-2*pi/3;pi/3];
    gimbal(:,i) = pt;
    d1(i) = point_to_line(pt',v11,v12);
    d2(i) = point_to_line(pt',v21,v22)  ;
    J = computeJac(pt,h0,beta);
    gDot(:,i) = J'/(J*J')*torque;
    Norm2(i) = norm(gDot(:,i),2);
    NormInf(i) = norm(gDot(:,i),"inf");
end

% Plots
figure
scatter3(gp(1,index1h),gp(2,index1h),gp(3,index1h), 10,'fill','b')
hold on
scatter3(gp(1,index3h),gp(2,index3h),gp(3,index3h), 10,'fill','g')
hold on
scatter3(gimbal(1,:),gimbal(2,:),gimbal(3,:), 10,'fill','r')
hold on
scatter3([-pi/3 0],[pi/3 -pi/3],[0 pi/3], 50,'fill','k')
hold on
plotcube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15)
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15)
zlabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15)
title('Gimbal trajectory','Interpreter','latex','FontSize',15)
legend('1h singularities','3h singularities','segment $r$','singular points','interpreter','latex','Location','southoutside','NumColumns',2)
grid off
axis equal
box on
text(-pi/3-0.5,pi/3+0.3,0.3,'A')
text(-0.5,-pi/3+0.3,pi/3+0.3,'B')
view(335,20)

figure
plot(indexvec,d1,'r-','LineWidth',1)
hold on
plot(indexvec,d2,'r-.','LineWidth',1)
xlabel('$k_r$','Interpreter','latex','FontSize',15)
ylabel('Euclidean distance [rad]','Interpreter','latex','FontSize',15)
title('Distance to the singularities','Interpreter','latex','FontSize',15)
legend('$d_A$','$d_B$','Interpreter','latex','FontSize',15)
grid off
box off

figure
% plot(indexvec,Norm2,'b-','LineWidth',1,'MarkerSize',2)
% hold on
plot(indexvec,NormInf,'b','LineWidth',1,'MarkerSize',2)
hold on
yline(1.5,'r','LineWidth',1)
xlabel('$k_r$','Interpreter','latex','FontSize',15)
ylabel('Norm [rad/s]','Interpreter','latex','FontSize',15);
xline(0,'k-.','LineWidth',0.5)
xline(1,'k-.','LineWidth',0.5)
% xline(-0.145,'m-.','LineWidth',1)
% xline(0.18,'m-.','LineWidth',1)
% xline(0.82,'m-.','LineWidth',1)
% xline(1.145,'m-.','LineWidth',1)
title('Largest control input','Interpreter','latex','FontSize',15)
legend('$|| \dot{\gamma} ||_{\infty}$',...
     '$\dot{\gamma}_{\max}$','Interpreter','latex','FontSize',15)
ylim([0 10])
grid off
box off

%% Functions

function h = MomentumCompute(g,h0,beta)
h = h0*[cos(beta + g(1)) + sin(g(2)) - cos(beta - g(3));
        sin(beta + g(1)) - cos(g(2)) + sin(beta - g(3))];
end

function J = computeJac(g,h0,beta)
J = h0*[-sin(beta + g(1))   cos(g(2))  -sin(beta - g(3));
         cos(beta + g(1))   sin(g(2))  -cos(beta - g(3))];
end

function d = point_to_line(pt, v1, v2)
      a = v2 - v1;
      b = pt - v1;
      d = norm(cross(a,b)) / norm(a);
end

