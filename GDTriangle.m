%%%%%%%%%%%%%%%%%%%%%%
% Triangle array
%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 30*pi/180;       % Triangle internal angle

%% Plot: Determinant
% Fix two axis and plot the determinant as function of the other two
n = 50;  % Number of points
g1Plot = linspace(-pi,pi,n);
g2Plot = linspace(-pi,pi,n);
g3Plot = 0;
D = []; G1 = []; G2 = [];
for i1 = 1:length(g1Plot)
    for i2 = 1:length(g2Plot)
        J = TriangleJacobian([g1Plot(i1) g2Plot(i2) g3Plot],h0,beta);
        D = [D det(J*J')];
        G1 = [G1 g1Plot(i1)];
        G2 = [G2 g2Plot(i2)];
    end
end
G1s = reshape(G1,n,n);
G2s = reshape(G2,n,n);
Ds = reshape(D,n,n);

figure
surf(Xs,Ys,Ds)
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\det$(JJ$^T$)','Interpreter','latex','FontSize',15);
title('Determinant visualization','Interpreter','latex','FontSize',15);
subtitle('$\gamma_3 = 0$','Interpreter','latex','FontSize',15);
xlim([-pi pi])
ylim([-pi pi])
grid off
axis square

%% Gradient
% Compute the first-order derivatives of the determinant
g = sym('g',[3 1]);
J = TriangleJacobian(g,h0,beta);
eq = det(J*J');
Dg1 = diff(eq,g(1));
Dg2 = diff(eq,g(2));
Dg3 = diff(eq,g(3));

%% Gradient descent
g1 = pi/3; g2 = 0; g3 = pi/4;   % Initial configuration
J = TriangleJacobian([g1 g2 g3],h0,beta);
disp(det(J*J'));

alpha = 0.2;        % Learning Rate
iterations = 20;    % Number of iterations
sigma = 0.01;       % Standard deviation (perturbations)
G1 = g1; G2 = g2; G3 = g3;
D = det(J*J');
for i = 1:iterations
    % Update 
    t1 = g1 - alpha*double(subs(Dg1));
    t2 = g2 - alpha*double(subs(Dg2));
    t3 = g3 - alpha*double(subs(Dg3));
    % Assign
    g1 = t1; g2 = t2; g3 = t3; 
    % Add perturbation
    g1 = g1 + normrnd(0,sigma);
    g2 = g2 + normrnd(0,sigma);
    g3 = g3 + normrnd(0,sigma);
    
    J = TriangleJacobian([g1 g2 g3],h0,beta);

    G1 = [G1 g1]; G2 = [G2 g2]; G3 = [G3 g3]; G4 = [G4 g4];
    D = [D det(J*J')];
end

%% Plot: Determinant Evolution
plot([0 1:iterations],D,'r','LineWidth',1)
xlabel('iterations','Interpreter','latex','FontSize',15);
ylabel('$\det$(JJ$^T$)','Interpreter','latex','FontSize',15);
title('Determinant evolution','Interpreter','latex','FontSize',15);
box off

%% Gimbal space
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

h1h = []; index1h = []; h3h = []; index3h = [];
for i = 1:length(gp)
    aux = TriangleMomentum(gp(:,i),h0,beta);
    if round(norm(aux)) == h0
        h1h = [h1h aux];
        index1h = [index1h i];
    else 
        h3h = [h3h aux];
        index3h = [index3h i];
    end
end

figure
scatter3(gp(1,index1h),gp(2,index1h),gp(3,index1h), 10,'fill','b')
hold on
scatter3(gp(1,index3h),gp(2,index3h),gp(3,index3h), 10,'fill','g')
hold on
plot3(G1,G2,G3,'r-*','LineWidth',1)
hold on
plotcube([2*pi 2*pi 2*pi],[-pi -pi -pi],0,[0 0 1]);
xlabel('$\gamma_1$ [rad]','Interpreter','latex','FontSize',15);
ylabel('$\gamma_2$ [rad]','Interpreter','latex','FontSize',15);
zlabel('$\gamma_3$ [rad]','Interpreter','latex','FontSize',15);
title('Singularities envelope','Interpreter','latex','FontSize',15);
subtitle('Triangular array','FontSize',12);
legend('1h singularities','3h singularities','Convergence','interpreter','latex','Location','southoutside','NumColumns',3);
axis equal
grid off
box on
view(335,20)

%% Function: Triangle Jacobian
function J = TriangleJacobian(g,h0,beta)
J = h0*[-sin(beta + g(1))   cos(g(2))  -sin(beta - g(3));
         cos(beta + g(1))   sin(g(2))  -cos(beta - g(3))];
end
%% Function: Triangle Momentum
function h = TriangleMomentum(g,h0,beta)
h = h0*[cos(beta + g(1)) + sin(g(2)) - cos(beta - g(3));
        sin(beta + g(1)) - cos(g(2)) + sin(beta - g(3))];
end