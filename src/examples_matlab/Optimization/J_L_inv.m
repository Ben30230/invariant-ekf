function [output] = J_L_inv(r)
%J_L_INV 此处显示有关此函数的摘要
%   此处显示详细说明
r=r(:);
if norm(r)<1e-11
    output = eye(3);
else
    output = eye(3)- 0.5*axis2skew(r)+(1/(r'*r)-(1+cos(norm(r)))/(2*norm(r)*sin(norm(r))))*axis2skew(r)*axis2skew(r);
end

end




