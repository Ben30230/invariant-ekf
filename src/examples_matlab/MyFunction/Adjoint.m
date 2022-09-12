function AdjX = Adjoint(X)
% Adjoint of SE_N(3)
N = size(X,2)-3;
R = X(1:3,1:3);
R_cell = repmat({R}, 1, N+1);
AdjX = blkdiag(R_cell{:});
for i = 1:N
    AdjX(3*i+1:3*i+3,1:3) = axis2skew(X(1:3,i+3))*R;
end
end