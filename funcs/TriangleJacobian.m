%% Function: Triangle Jacobian
function J = TriangleJacobian(g,h0,lambda)
J = h0*[-sin(lambda + g(1))   cos(g(2))  -sin(lambda - g(3));
         cos(lambda + g(1))   sin(g(2))  -cos(lambda - g(3))];
end