addpath('functionsformytest\')
addpath('MyFunction\')
addpath('Optimization\')
%%
clc,clear
w_true = rand(3,1000);
perturb = 0.001*rand(3,1);
Delta_t = 0.001;
R0 = eye(3);

R_true = R0;
R_seprate = R0;
for i=1:length(w_true)
    R_true = R_true * Exp((w_true(:,i)-perturb)*Delta_t);
    R_seprate = R_seprate * Exp(w_true(:,i)*Delta_t)*Exp(-perturb*Delta_t);
end
R_true
R_seprate
norm(R_true'*R_seprate)


