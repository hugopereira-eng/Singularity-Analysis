%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo simulation
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
D = simplify(det(J*J'));
G = gradient(D);
H = hessian(D);

%% Random set
% Uniformly generated samples
samples = 5;
S = zeros(4,samples);
ind = 0;
for i = 1:samples
    S(:,i) = [-pi + 2*pi*rand;
              -pi + 2*pi*rand;
              -pi + 2*pi*rand;
              -pi + 2*pi*rand];
end

%% Newton
I = zeros(1,samples);
hbar = waitbar(0,'Simulation Progress');
for i = 1:length(S)
    g = S(:,i);
    g1 = g(1); g2 = g(2); g3 = g(3); g4 = g(4);
    J = PyramidJacobian(g,h0,beta);
    D = det(J*J');
    iterations = 0;
    while D(iterations+1) > 10e-3
        g = [g1;g2;g3;g4];
        % Update
        g = g - (double(subs(H))+diag([0.1 0.1 0.1 0.1]))\double(subs(G));
       
        J = PyramidJacobian(g,h0,beta);
        D = [D det(J*J')];
        
        g1 = g(1); g2 = g(2); g3 = g(3); g4 = g(4);
        iterations = iterations + 1;
    end
    I(i) = iterations;
    waitbar(i/length(S), hbar);
    disp(i)
end
close(hbar)

%% Results
plot(1:samples,I,'r','LineWidth',1)
xlabel('samples','Interpreter','latex','FontSize',15);
ylabel('iterations','Interpreter','latex','FontSize',15);
title('Iterations','Interpreter','latex','FontSize',15);
box off
fprintf("Mean: %d \n",mean(I));
fprintf("Max: %d \n",max(I));
