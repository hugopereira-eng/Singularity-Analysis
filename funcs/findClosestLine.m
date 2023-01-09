%% Function: Find the closest singularity line
function p = FindClosestLine(g)
% list   = [   -pi,  2*pi/3, -2*pi/3;
%              -pi,  2*pi/3,    pi/3;
%              -pi,  2*pi/3, -5*pi/3;
%              -pi,  2*pi/3,  4*pi/3;
%              -pi,   -pi/3,    pi/3;
%              -pi,   -pi/3, -2*pi/3;
%            -pi/3,    pi/3,     -pi;
%          -2*pi/3,     -pi,  2*pi/3;
%          -2*pi/3,     -pi,  - pi/3;
%            -pi/3, -2*pi/3,     -pi;
%           2*pi/3,    pi/3,     -pi;
%             pi/3,     -pi,  2*pi/3;
%             pi/3,     -pi,   -pi/3;
%           2*pi/3, -2*pi/3,     -pi;];

list   = [   -pi/3,     pi/3,  0;
             -pi/3,   4*pi/3,  0;
             -pi/3,   7*pi/3,  0;
             -pi/3,  -2*pi/3,  0;
             -pi/3,  -5*pi/3,  0;

            2*pi/3,     pi/3,  0;
            2*pi/3,   4*pi/3,  0;
            2*pi/3,   7*pi/3,  0;
            2*pi/3,  -2*pi/3,  0;
            2*pi/3,  -5*pi/3,  0;
             
            5*pi/3,     pi/3,  0;
            5*pi/3,   4*pi/3,  0;
            5*pi/3,   7*pi/3,  0;
            5*pi/3,  -2*pi/3,  0;
            5*pi/3,  -5*pi/3,  0;
             
           -4*pi/3,     pi/3,  0;
           -4*pi/3,   4*pi/3,  0;
           -4*pi/3,   7*pi/3,  0;
           -4*pi/3,  -2*pi/3,  0;
           -4*pi/3,  -5*pi/3,  0;
             
           -7*pi/3,     pi/3,  0;
           -7*pi/3,   4*pi/3,  0;
           -7*pi/3,   7*pi/3,  0;
           -7*pi/3,  -2*pi/3,  0;
           -7*pi/3,  -5*pi/3,  0;];

v = [1;1;1];   % Vector director
dmin = 9999;
for i = 1:length(list)
    d = PointToLine(g, list(i,:)', list(i,:)' + v);
    if (d < dmin)
        dmin = d;
        p = list(i,:)';
    end
end
end

function d = PointToLine(pt, v1, v2)
      a = v2 - v1;
      b = pt - v1;
      d = norm(cross(a,b)) / norm(a);
end