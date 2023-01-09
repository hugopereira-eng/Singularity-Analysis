%% Function: Triangle momentum
function h = TriangleMomentum(g,h0,lambda)
    h = h0*[cos(lambda+g(1)) + sin(g(2)) - cos(lambda-g(3));
            sin(lambda+g(1)) - cos(g(2)) + sin(lambda-g(3))];
end