%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-CMG determinant expression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG

%% Symbolic expression
g1 = sym('g1', 'real');
g2 = sym('g2', 'real');
g3 = sym('g3', 'real');
g4 = sym('g4', 'real');
% beta = sym('beta', 'real');
g = [g1,g2,g3,g4];  
J = PyramidJacobian(g,h0,beta);
% J = RoofJacobian(g,h0);
disp(simplify(det(J*J')))
% r = [J(:,1) J(:,2) J(:,3)];
% disp(simplify(det(r)^2))

%% Gradient
D = simplify(det(J*J'));
disp(D);
G = gradient(D);
H = hessian(D);

%% Analytical solution
D1 = sin(beta)*(sin(g3)*sin(g4)*cos(g1) + cos(g3)*sin(g4)*sin(g1) + ...
    cos(beta)*(cos(g3)*cos(g4)*sin(g1) - sin(g3)*cos(g4)*cos(g1)) + ...
    2*cos(beta)^2*cos(g3)*cos(g4)*cos(g1));

D2 = sin(beta)*(sin(g2)*sin(g3)*cos(g4) + cos(g2)*sin(g3)*sin(g4) + ...
    cos(beta)*(cos(g2)*cos(g3)*sin(g4) - sin(g2)*cos(g3)*cos(g4)) + ...
    2*cos(beta)^2*cos(g2)*cos(g3)*cos(g4));
% D1 = det(J(1:2,1:2));
% D1 = det([J(1:2,1) J(1:2,4)]);
% D2 = det([J(1:2,2) J(1:2,4)]);
% D3 = det(J(1:2,3:4));

eq1 = D1 == 0;
eq2 = D2 == 0;
% eq3 = D3 == 0;
% eq4 = D4 == 0;
% S = solve(eq1,eq2,g1,g2,"Real",true,"PrincipalValue",false,"ReturnConditions",true);
S = solve(eq1,eq2,g1,g2,"Real",true,"PrincipalValue",false);
% disp(vpa(simplify(S.g1),5))
% disp(vpa(simplify(S.g2),5))

%% Parameters
beta = 54.73*pi/180;    % Pyramid skew angle

%% Test 
g = [90; 45; 45; 90]*pi/180;
J = PyramidJacobian(g,h0,beta);
% J = RoofJacobian(g,h0);
fprintf("Determinant: %d \n",det(J*J'));

%% Gamma - Solve equations
g1 = sym('g1', 'real');
g2 = sym('g2', 'real');

n = 30;
g3 = linspace(-pi,pi,n);
g4 = linspace(-pi,pi,n);
% g4 = 180*pi/180;
s1 = []; s2 = []; g3Plot = []; g4Plot = [];
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
        S = solve(eq1,eq2,g1,g2);
        s1 = [s1 real(double(S.g1))];
        s2 = [s2 real(double(S.g2))];
        g3Plot = [g3Plot g3(i3)];
        g4Plot = [g4Plot g4(i4)];
    end
    disp(i3)
end

%% Gamma - Use solutions
n3 = 100;
n4 = 1;
g3v = linspace(-pi,pi,n3);
% g4v = linspace(-pi,pi,n4);
g4v = 0*pi/180;
s1 = []; s2 = []; 
g3Plot = []; g4Plot = [];
t1 = S.g1(1);
t2 = S.g2(1);
for j1 = 1:length(S.g1)
    t1 = S.g1(j1);
    for j2 = 1:length(S.g2)
        t2 = S.g2(j2);
        for i4 = 1:length(g4v)
            g4 = g4v(i4);
            for i3 = 1:length(g3v)
                g3 = g3v(i3);
                g1 = subs(t1);
                g2 = subs(t2);
                s1 = [s1 real(double(g1))];
                s2 = [s2 real(double(g2))];
                g3Plot = [g3Plot g3v(i3)];
                g4Plot = [g4Plot g4v(i4)];
            end
        end
        disp(j2)
    end
end

%% Plot - Gamma plot (g1)
figure
g1Plot = s1;
C = sqrt(g1Plot.^2+g3Plot.^2+g4Plot.^2);
scatter3(g1Plot,g3Plot,g4Plot,10,C,'fill')
hold on
PlotCube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\gamma_4$ [rad]','Interpreter','latex','FontSize',15);
title('Singularity envelope','Interpreter','latex','FontSize',15);
subtitle('Without CMG \#2','Interpreter','latex','FontSize',12);
xlim([-pi pi])
ylim([-pi pi])
zlim([-pi pi])
grid off
axis equal
box on
% view(155,20)
view(335,20)

%% Plot - Gamma plot (g2)
figure
g2Plot = s2;
C = sqrt(g2Plot.^2+g3Plot.^2+g4Plot.^2);
scatter3(g2Plot,g3Plot,g4Plot,10,C,'fill')
hold on
PlotCube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\gamma_4$ [rad]','Interpreter','latex','FontSize',15);
title('Singularity envelope','Interpreter','latex','FontSize',15);
subtitle('Without CMG \#1','Interpreter','latex','FontSize',12);
xlim([-pi pi])
ylim([-pi pi])
zlim([-pi pi])
grid off
axis equal
box on
% view(155,20)
view(335,20)

%% Plot - Gamma plot (g4 fixed)
figure
g1Plot = s1;
g2Plot = s2;
scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
hold on
PlotCube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
title('Singularity envelope','Interpreter','latex','FontSize',15);
subtitle('$\gamma_4 = 0^{\circ}$','Interpreter','latex','FontSize',12);
xlim([-pi pi])
ylim([-pi pi])
zlim([-pi pi])
grid off
axis equal
box on
view(335,20)

%% Plot - Gamma plot (colormap) as a function of g4
figure
for j = n4/2:n4
    g1Plot = []; g2Plot = []; g3Plot_ = [];
    for i = 1:16
        g1Plot  = [g1Plot  s1((n3*(j-1)+1+n3*n4*(i-1)):(n3*j+n3*n4*(i-1)))];
        g2Plot  = [g2Plot  s2((n3*(j-1)+1+n3*n4*(i-1)):(n3*j+n3*n4*(i-1)))];
        g3Plot_ = [g3Plot_ g3Plot((n3*(j-1)+1+n3*n4*(i-1)):(n3*j+n3*n4*(i-1)))]; 
    end
    C = ones(16*n3,1)*g4v(j);
    scatter3(g1Plot,g2Plot,g3Plot_,10,C,'fill')
    hold on
end
PlotCube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
title('Singularity envelope','Interpreter','latex','FontSize',15);
subtitle('$\gamma_4 \in [0,\pi]$','Interpreter','latex','FontSize',12);
xlim([-pi pi])
ylim([-pi pi])
zlim([-pi pi])
grid off
axis equal
box on
colorbar
view(335,20)


% %% Plot - Gamma (extra)
% figure
% g1Plot = s1(1,:);
% g2Plot = s2(1,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(1,:);
% g2Plot = s2(2,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(1,:);
% g2Plot = s2(3,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(1,:);
% g2Plot = s2(4,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(2,:);
% g2Plot = s2(1,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(2,:);
% g2Plot = s2(2,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(2,:);
% g2Plot = s2(3,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(2,:);
% g2Plot = s2(4,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(3,:);
% g2Plot = s2(1,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(3,:);
% g2Plot = s2(2,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(3,:);
% g2Plot = s2(3,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(3,:);
% g2Plot = s2(4,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(4,:);
% g2Plot = s2(1,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(4,:);
% g2Plot = s2(2,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(4,:);
% g2Plot = s2(3,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% g1Plot = s1(4,:);
% g2Plot = s2(4,:);
% scatter3(g1Plot,g2Plot,g3Plot, 10,'fill','b')
% hold on
% PlotCube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
% xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
% ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
% zlabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
% title('Singularity envelope','Interpreter','latex','FontSize',15);
% xlim([-pi pi])
% ylim([-pi pi])
% zlim([-pi pi])
% grid off
% axis equal
% box on
% 
