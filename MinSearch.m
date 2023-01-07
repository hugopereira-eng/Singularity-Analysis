
h0 = 1;
beta = 54.73*pi/180;

g0 = [-pi -pi 0 0];
tic
sol = fminsearch(@fun,g0);
toc
J = PyramidJacobian(sol,h0,beta);
disp(det(J*J'))


%% Det

function D = fun(g)
    h0 = 1;
    beta = 54.73*pi/180;
    J = PyramidJacobian(g,h0,beta);
    D = det(J*J');
end

%% Jacobian
function J = PyramidJacobian(g,h0,beta)

J = h0*[-cos(beta)*cos(g(1)) sin(g(2)) cos(beta)*cos(g(3)) -sin(g(4));
        -sin(g(1)) -cos(beta)*cos(g(2)) sin(g(3)) cos(beta)*cos(g(4));
        sin(beta)*cos(g(1)) sin(beta)*cos(g(2)) sin(beta)*cos(g(3)) sin(beta)*cos(g(4))];
end
