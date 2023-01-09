%% Function: Distance between a point and a line
function d = Point2Line(pt, v1, v2)
      a = v2 - v1;
      b = pt - v1;
      d = norm(cross(a,b)) / norm(a);
end