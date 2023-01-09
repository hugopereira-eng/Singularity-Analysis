function J = TriangleJacobian(g,h0,alpha)
J = h0*[-sin(alpha + g(1))   cos(g(2))  -sin(alpha - g(3));
         cos(alpha + g(1))   sin(g(2))  -cos(alpha - g(3))];
end