%% Function: Roof momentum
function h = RoofMomentum(g,h0)
    h = h0*[ sin(g(1)) + cos(g(2));
            -cos(g(3)) + sin(g(4));
            -cos(g(1)) + sin(g(2)) + sin(g(3)) + cos(g(4))];
end