%% Function: Pyramid momentum
function h = PyramidMomentum(g,h0,beta)
    h = h0*[-cos(beta)*sin(g(1)) - cos(g(2))           + cos(beta)*sin(g(3)) + cos(g(4));
             cos(g(1))           - cos(beta)*sin(g(2)) - cos(g(3))           + cos(beta)*sin(g(4));
             sin(beta)*sin(g(1)) + sin(beta)*sin(g(2)) + sin(beta)*sin(g(3)) + sin(beta)*sin(g(4))];
end