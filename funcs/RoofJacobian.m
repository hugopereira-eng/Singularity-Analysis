%% Function: Roof Jacobian
function J = RoofJacobian(g,h0)
J = h0*[ cos(g(1))   -sin(g(2))             0            0;
                 0            0     sin(g(3))    cos(g(4));
         sin(g(1))    cos(g(2))     cos(g(3))   -sin(g(4))];
end