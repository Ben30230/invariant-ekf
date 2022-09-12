function [output] = Log(X)
%LOG 此处显示有关此函数的摘要
%   此处显示详细说明
% output=skew2axis(logm(X));
theta = acos((trace(X)-1)/2);
if abs(theta) < 1e-11
    output = zeros(3,1);
else
    w = [X(3,2)-X(2,3);X(1,3)-X(3,1);X(2,1)-X(1,2)]/(2*sin(theta));
    output = w * theta;
end



end

