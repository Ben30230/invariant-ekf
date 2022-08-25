function [output] = J_R(r)
%J_R r \in R^n
%
if norm(r)<1e-11
    output = eye(3);
else
    output = eye(3)-(1-cos(norm(r)))/(r'*r)*axis2skew(r)+(norm(r)-sin(norm(r)))/(norm(r)^3)*(axis2skew(r)^2) ;
end


end

