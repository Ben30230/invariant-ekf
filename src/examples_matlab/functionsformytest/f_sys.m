function [y] = f_sys(x)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
y=zeros(3,1);

y(1) = x(1)+0.5*x(2) +x(1)*x(3);
y(2) = 2*x(1)+x(2)+x(2)*x(3);
y(3) = x(3);


end