function [p, q, p_cross_q] = plucker_line(x, Q, s)
p = Q*x;
p = p/norm(p(:));
q = cross(s, p);
p_cross_q = cross(p, q);