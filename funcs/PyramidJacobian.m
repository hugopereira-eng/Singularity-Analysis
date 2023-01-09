%% Function: Pyramid Jacobian
function J = PyramidJacobian(g,h0,beta)
J = h0*[-cos(beta)*cos(g(1))   sin(g(2))             cos(beta)*cos(g(3))  -sin(g(4));
        -sin(g(1))            -cos(beta)*cos(g(2))   sin(g(3))             cos(beta)*cos(g(4));
         sin(beta)*cos(g(1))   sin(beta)*cos(g(2))   sin(beta)*cos(g(3))   sin(beta)*cos(g(4))];
end  