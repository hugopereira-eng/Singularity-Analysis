close all
clear
clc

%% Parameters

h0 = 1;

%% Symbolic determinant expression

g1 = sym('g1', 'real');
g2 = sym('g2', 'real');
g3 = sym('g3', 'real');
g4 = sym('g4', 'real');
g = [g1,g2,g3,g4];  
J = computeJac(g,h0);
disp(simplify(det(J*J')))
r = [J(:,1) J(:,2) J(:,3)];
disp(simplify(det(r)^2))

%% Try values

g = [0; 45; 180; 0]*pi/180;
J = computeJac(g,h0);
fprintf("Determinant: %d \n",det(J*J'));

%% Functions

function J = computeJac(g,h0)
J = h0*[ cos(g(1))   -sin(g(2))             0            0;
                 0            0     sin(g(3))    cos(g(4));
         sin(g(1))    cos(g(2))     cos(g(3))   -sin(g(4))];
end