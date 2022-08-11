function [y] = multirandn(mu,Sigma)
% generate y ~ N(mu,Sigma)
%   If x ~ N(0,I)    we have y=cx + mu ~ N(mu,cIc')

mu=mu(:);
n=length(mu);

x=randn(n,1);
c=Sigma^(1/2);

y=c*x+mu;



end