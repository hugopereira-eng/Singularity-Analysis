%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-CMG determinant expression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("..\funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 54.73*pi/180;    % Pyramid skew angle

%% Symbolic expression
g1 = sym('g1', 'real');
g2 = sym('g2', 'real');
g3 = sym('g3', 'real');
g4 = sym('g4', 'real');
g = [g1,g2,g3,g4];  
J = PyramidJacobian(g,h0,beta);
% J = RoofJacobian(g,h0);
disp(simplify(det(J*J')))
r = [J(:,1) J(:,2) J(:,3)];
disp(simplify(det(r)^2))

%% Test 
g = [90; 45; 45; 90]*pi/180;
J = PyramidJacobian(g,h0,beta);
% J = RoofJacobian(g,h0);
fprintf("Determinant: %d \n",det(J*J'));
