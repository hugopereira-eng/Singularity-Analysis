%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hugo Pereira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("funcs\")

%% Parameters
h0 = 1;                  % Angular momentum of each CMG
beta = 30*pi/180;        % Triangle inner angle

%% Gradient
g = sym('g',[3 1]);
J = TriangleJacobian(g,h0,beta);
eq = det(J*J');
Dg1 = diff(eq,g(1));
Dg2 = diff(eq,g(2));
Dg3 = diff(eq,g(3));

%% Random set
% Uniformly generated samples
n = 10;
samples = n^3;
G = zeros(3,samples);
ind = 0;
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            ind = ind + 1;
            G(:,ind) = [-pi/2 + pi*rand;
                        -pi/2 + pi*rand;
                        -pi/2 + pi*rand];
        end
    end
end

%% Gradient descent
alpha = 0.05;         % learning rate
sigma = 0.0001;       % Standard deviation (perturbations)
d = zeros(1,samples);
GD = zeros(3,samples);
hbar = waitbar(0,'Simulation Progress');
for i = 1:length(G)
    g1 = G(1,i); g2 = G(2,i); g3 = G(3,i);
    J = TriangleJacobian([g1 g2 g3],h0,beta);
    D = det(J*J');
    iterations = 0;
    while D(iterations+1) > 10e-5
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
        D = [D det(J*J')];

        iterations = iterations + 1;
    end
    GD(:,i) = [g1;g2;g3];
    waitbar(i/length(G), hbar);
end
close(hbar)

%% Closest

v = [1;1;1];                         % Vector director
As  = v*v'/(v'*v);                   % Projection matrix
GC = zeros(3,samples);
for i = 1:length(G)
    g0 = findClosestLine(G(:,i));
    bs  = (eye(3) -As)*g0;        % Offset vector 
    gs0 = As*G(:,i) + bs;
    GC(:,i) = gs0;
end

%% Error

error = GD - GC;
e = zeros(1,samples);
for i = 1:length(G)
    e(i) = norm(error(:,i),2);
end
plot(e)
