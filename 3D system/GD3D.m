%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gradient Descent with the 3D system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 54.73*pi/180;    % Pyramid skew angle

%% Plot: Determinant manifold
% Fix two axis and plot the determinant as function of the other two
n = 50;  % Number of points
g1Plot = linspace(-pi,pi,n);
g2Plot = linspace(-pi,pi,n);
g3Plot = 0;
g4Plot = 0;
D = []; G1 = []; G2 = [];
for i1 = 1:length(g1Plot)
    for i2 = 1:length(g2Plot)
        J = PyramidJacobian([g1Plot(i1) g2Plot(i2) g3Plot g4Plot],h0,beta);
%         J = RoofJacobian([g1Plot(i1) g2Plot(i2) g3Plot g4Plot],h0);
        D = [D det(J*J')];
        G1 = [G1 g1Plot(i1)];
        G2 = [G2 g2Plot(i2)];
    end
end
G1s = reshape(G1,n,n);
G2s = reshape(G2,n,n);
Ds = reshape(D,n,n);

figure
surf(G1s,G2s,Ds)
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\det$(JJ$^T$)','Interpreter','latex','FontSize',15);
title('Determinant visualization','Interpreter','latex','FontSize',15);
subtitle('$\gamma_3 = \gamma_4 = 0$','Interpreter','latex','FontSize',15);
xlim([-pi pi])
ylim([-pi pi])
grid off
axis square

%% Gradient
% Compute the first-order derivatives of the determinant
g = sym('g',[4 1]);
J = PyramidJacobian(g,h0,beta);
% J = RoofJacobian(g,h0);
eq = det(J*J');
Dg1 = diff(eq,g(1));
Dg2 = diff(eq,g(2));
Dg3 = diff(eq,g(3));
Dg4 = diff(eq,g(4));

%% Initial configuration
g1 = 0; g2 = 0; g3 = pi/2; g4 = pi/4;
J = PyramidJacobian([g1 g2 g3 g4],h0,beta);
% J = RoofJacobian([g1 g2 g3 g4],h0);
disp(det(J*J'));

%% Gradient descent
alpha = 0.1;         % Learning rate
iterations = 20;     % Number of iterations
sigma = 0.001;       % Standard deviation (perturbations)
G1 = g1; G2 = g2; G3 = g3; G4 = g4; 
D = det(J*J');
for i = 1:iterations
    % Update 
    t1 = g1 - alpha*double(subs(Dg1));
    t2 = g2 - alpha*double(subs(Dg2));
    t3 = g3 - alpha*double(subs(Dg3));
    t4 = g4 - alpha*double(subs(Dg4));
    % Assign
    g1 = t1; g2 = t2; g3 = t3; g4 = t4;
    % Add perturbation
    g1 = g1 + normrnd(0,sigma);
    g2 = g2 + normrnd(0,sigma);
    g3 = g3 + normrnd(0,sigma);
    g4 = g4 + normrnd(0,sigma);
    
    J = PyramidJacobian([g1 g2 g3 g4],h0,beta);
%     J = RoofJacobian([g1 g2 g3 g4],h0);

    G1 = [G1 g1]; G2 = [G2 g2]; G3 = [G3 g3]; G4 = [G4 g4];
    D = [D det(J*J')];
end

%% Plot: Determinant Evolution
plot([0 1:iterations],D,'r','LineWidth',1)
xlabel('iterations','Interpreter','latex','FontSize',15);
ylabel('$\det$(JJ$^T$)','Interpreter','latex','FontSize',15);
title('Determinant evolution','Interpreter','latex','FontSize',15);
box off
