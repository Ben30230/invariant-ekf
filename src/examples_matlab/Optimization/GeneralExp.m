function [output] = GeneralExp(xi)
%GENERALEXP 此处显示有关此函数的摘要
%   此处显示详细说明
xi=xi(:);
n = length(xi)/3;
output = eye(n+2);
output(1:3,1:3) = Exp(xi(1:3));
for i = 4:n+2
    output(1:3,i) = J_L(xi(1:3)) * xi((i-2)*3-2:(i-2)*3);
end


end

