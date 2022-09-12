function [X_inv] = Liegroup_inverse(X)
%LIEGROUP_INVERSE 此处显示有关此函数的摘要
%   此处显示详细说明
n = size(X,2);
X_inv = eye(size(X));
R_T = X(1:3,1:3)';
X_inv(1:3,1:3) = R_T;
for i = 4:n
    X_inv(1:3,i) = -R_T * X(1:3,i);
end

end

