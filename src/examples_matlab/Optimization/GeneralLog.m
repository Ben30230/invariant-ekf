function [output] = GeneralLog(X)
%GENERALLOG 此处显示有关此函数的摘要
%   此处显示详细说明
n = size(X,2)-2;
output = zeros(n*3,1);
phi = Log(X(1:3,1:3));
output(1:3) = phi;

for i=2:n
    output(i*3-2:i*3) = J_L_inv(phi) * X(1:3,i+2);
end


end

