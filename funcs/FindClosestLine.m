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
    d = Point2Line(g, list(i,:)', list(i,:)' + v);
    if (d < dmin)
        dmin = d;
        p = list(i,:)';
    end
end
end