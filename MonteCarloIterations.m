%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("funcs\")

%% Parameters
h0 = 1;                 % Angular momentum of each CMG
beta = 54.73*pi/180;    % Pyramid skew angle

%% Gradient & Hessian
g = sym('g',[4 1]);
J = PyramidJacobian(g,h0,beta);
eq = det(J*J');
Dg = gradient(eq);
Hg = hessian(eq);

%% Random set
% Uniformly generated samples
samples = 10;
G = zeros(4,samples);
ind = 0;
for i = 1:samples
    G(:,i) = [-pi + 2*pi*rand;
              -pi + 2*pi*rand;
              -pi + 2*pi*rand;
              -pi + 2*pi*rand];
end

%% Gradient descent
sigma = 0.001;       % Standard deviation (perturbations)
I = zeros(1,samples);
hbar = waitbar(0,'Simulation Progress');
for i = 1:length(G)
    g1 = G(1,i); g2 = G(2,i); g3 = G(3,i); g4 = G(4,i);
    g = [g1;g2;g3;g4];
    J = PyramidJacobian(g,h0,beta);
    D = det(J*J');
    iterations = 0;
    while D(iterations+1) > 10e-3

        % alpha
        Dv = double(subs(Dg));
        Hv = double(subs(Hg));
        gradD = Dv - Hv*g;
        alpha = gradD'*gradD/(gradD'*Hv*gradD);

        % Update
        g = g - alpha*gradD;

%         % Update
%         t1 = g1 - alpha*double(subs(Dg1));
%         t2 = g2 - alpha*double(subs(Dg2));
%         t3 = g3 - alpha*double(subs(Dg3));
%         t4 = g4 - alpha*double(subs(Dg4));
%         % Add perturbation - generated from a Gaussian distribution
%         g1 = t1 + normrnd(0,sigma); 
%         g2 = t2 + normrnd(0,sigma);
%         g3 = t3 + normrnd(0,sigma); 
%         g4 = t4 + normrnd(0,sigma);
        
        J = PyramidJacobian(g,h0,beta);
        D = [D det(J*J')];

        iterations = iterations + 1;
    end
    I(i) = iterations;
    waitbar(i/length(G), hbar);
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
