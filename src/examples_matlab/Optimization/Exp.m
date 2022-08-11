function [output] = Exp(xi)
%EXP xi \in R^{3n}
% 
% xi_l=length(xi)/3;
% if mod(xi_l,1)~=0
%     error("dim of input error!");
% end


output=expm(axis2skew(xi));

end

