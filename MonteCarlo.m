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
eq = det(J*J');
Dg1 = diff(eq,g(1));
Dg2 = diff(eq,g(2));
Dg3 = diff(eq,g(3));
Dg4 = diff(eq,g(4));

%% Random set
% Uniformly generated samples
n = 5;
samples = n^4;
G = zeros(4,samples);
ind = 0;
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            for i4 = 1:n
                ind = ind + 1;
                G(:,ind) = [-pi + 2*pi*rand;
                            -pi + 2*pi*rand;
                            -pi + 2*pi*rand;
                            -pi + 2*pi*rand];
            end
        end
    end
end

%% Gradient descent
alpha = 0.1;         % learning rate
sigma = 0.001;       % Standard deviation (perturbations)
I = zeros(1,samples);
hbar = waitbar(0,'Simulation Progress');
for i = 1:length(G)
    g1 = G(1,i); g2 = G(2,i); g3 = G(3,i); g4 = G(4,i);
    J = PyramidJacobian([g1 g2 g3 g4],h0,beta);
    D = det(J*J');
    iterations = 0;
    while D(iterations+1) > 10e-3
        % Update
        t1 = g1 - alpha*double(subs(Dg1));
        t2 = g2 - alpha*double(subs(Dg2));
        t3 = g3 - alpha*double(subs(Dg3));
        t4 = g4 - alpha*double(subs(Dg4));
        % Add perturbation - generated from a Gaussian distribution
        g1 = t1 + normrnd(0,sigma); 
        g2 = t2 + normrnd(0,sigma);
        g3 = t3 + normrnd(0,sigma); 
        g4 = t4 + normrnd(0,sigma);
        
        J = PyramidJacobian([g1 g2 g3 g4],h0,beta);
        D = [D det(J*J')];

        iterations = iterations + 1;
    end
    I(i) = iterations;
    waitbar(i/length(G), hbar);
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
