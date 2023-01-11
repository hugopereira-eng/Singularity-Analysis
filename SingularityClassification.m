%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Singularity classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 54.73*pi/180;    % Pyramid skew angle

%% Gradient
g = sym('g',[4 1]);
J = PyramidJacobian(g,h0,beta);
eq = det(J*J');
Dg1 = diff(eq,g(1));
Dg2 = diff(eq,g(2));
Dg3 = diff(eq,g(3));
Dg4 = diff(eq,g(4));

G = [Dg1; Dg2; Dg3; Dg4];

%% Hessian
Dg11 = diff(Dg1,g(1));
Dg12 = diff(Dg1,g(2));
Dg13 = diff(Dg1,g(3));
Dg14 = diff(Dg1,g(4));

Dg21 = diff(Dg2,g(1));
Dg22 = diff(Dg2,g(2));
Dg23 = diff(Dg2,g(3));
Dg24 = diff(Dg2,g(4));

Dg31 = diff(Dg3,g(1));
Dg32 = diff(Dg3,g(2));
Dg33 = diff(Dg3,g(3));
Dg34 = diff(Dg3,g(4));

Dg41 = diff(Dg4,g(1));
Dg42 = diff(Dg4,g(2));
Dg43 = diff(Dg4,g(3));
Dg44 = diff(Dg4,g(4));

H = [Dg11 Dg12 Dg13 Dg14;
     Dg21 Dg22 Dg23 Dg24;
     Dg31 Dg32 Dg33 Dg34;
     Dg41 Dg42 Dg43 Dg44;];


%% Evaluate
g1 = 0; g2 = -2*pi/3; g3 = pi/6; g4 = -2*pi/3; 

Gv = double(subs(G));
Hv = double(subs(H));

disp(Gv)
disp(Hv)
disp(eig(Hv))
