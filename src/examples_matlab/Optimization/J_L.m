function [output] = J_L(r)
%J_L 此处显示有关此函数的摘要
%   此处显示详细说明

if norm(r)<1e-11
    output = eye(3);
else
    output = eye(3)+(1-cos(norm(r)))/(r'*r)*axis2skew(r)+(norm(r)-sin(norm(r)))/(norm(r)^3)*(axis2skew(r)^2) ;
end

end

