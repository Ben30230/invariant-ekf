function [output] = Exp(xi)
%EXP xi \in R^{3n}
%
% xi_l=length(xi)/3;
% if mod(xi_l,1)~=0
%     error("dim of input error!");
% end


% output=expm(axis2skew(xi));

% exp_w(xi/norm(xi),norm(xi));

if sum(abs(xi))<1e-15
    output = eye(3);
else
    output = eye(3) + axis2skew(xi)/norm(xi)*sin(norm(xi))+axis2skew(xi)^2/(xi'*xi)*(1-cos(norm(xi)));
end

end

