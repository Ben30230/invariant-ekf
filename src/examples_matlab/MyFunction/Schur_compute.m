function [sol,A_schur,b_schur] = Schur_compute(H,z)
% Hx = z solve for x(half:end)
%

% output --- A_schur *sol = b_schur

% if ~issymmetric(H)
%     error("not a symmetric matrix")
% end

n = size(H,2);

A = H(1:n/2,1:n/2);
B = H(1:n/2,n/2+1:end);
BT = H(n/2+1:end,1:n/2);
C = H(n/2+1:end,n/2+1:end);

A_schur = C-BT/A*B;
b_schur = z(n/2+1:end)-BT/A*z(1:n/2);
sol = A_schur\b_schur;
end

