function p = findClosestLine(g)

list   = [   -pi,  2*pi/3, -2*pi/3;
             -pi,  2*pi/3,    pi/3;
             -pi,   -pi/3,    pi/3;
             -pi,   -pi/3, -2*pi/3;
           -pi/3,    pi/3,     -pi;
         -2*pi/3,     -pi,  2*pi/3;
         -2*pi/3,     -pi,  - pi/3;
           -pi/3, -2*pi/3,     -pi;
          2*pi/3,    pi/3,     -pi;
            pi/3,     -pi,  2*pi/3;
            pi/3,     -pi,   -pi/3;
          2*pi/3, -2*pi/3,     -pi;];

v = [1;1;1];   % vector director of the singularities
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