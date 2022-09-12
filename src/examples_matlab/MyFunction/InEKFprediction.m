function [output] = InEKFprediction(Xk,Delta_t,w_k,acc_k,bg_bar,ba_bar)
%INEKFPREDICTION 此处显示有关此函数的摘要
%   此处显示详细说明
gravity = [0 0 -9.81]';

if nargin < 5 
    bg_bar = zeros(3,1);
    ba_bar = zeros(3,1);
end

Gamma_i = [eye(3),Delta_t*gravity,0.5*gravity*Delta_t*Delta_t,zeros(3,1);zeros(3),eye(3)];
Phi_i = Xk;
Phi_i(1:3,5) = Phi_i(1:3,5)+Delta_t * Phi_i(1:3,4);
Upsilon_i = [Exp((w_k-bg_bar)*Delta_t),(acc_k-ba_bar)*Delta_t,zeros(3,2);zeros(3),eye(3)];
output = Gamma_i*Phi_i*Upsilon_i;

end

