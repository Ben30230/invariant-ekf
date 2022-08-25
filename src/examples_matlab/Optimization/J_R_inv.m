function [output] = J_R_inv(r)
%J_R_INV r \in R^n
%   
if norm(r)<1e-11
    output = eye(3);
else
    output = eye(3) + 0.5*axis2skew(r)+(1/norm(r)^2-(1+cos(norm(r)))/(2*norm(r)*sin(norm(r))))*(axis2skew(r)^2);
end


end

