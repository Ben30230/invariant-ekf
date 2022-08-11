%% Example 1
clc,clear
% target
f_T = [0,0,0]';
R_hat = angvec2r(pi/2,[0,0,1]);

% function
f=@(r)Log(R_hat'*Exp(r)); 

% jacobian 
J_Rk=@(r)(eye(3)-(1-cos(norm(r)))/(r'*r)*axis2skew(r)+(norm(r)-sin(norm(r)))/(norm(r)^3)*(axis2skew(r)^2)) ;
J_k = @(r)eye(3);

R = Exp(rand(3,1));

for i =1:100
    tau = Log(R);
    epsilon = f(tau)-f_T;
    delta_R = -((J_k(tau)/J_Rk(tau))'*(J_k(tau)/J_Rk(tau)))\((J_k(tau)/J_Rk(tau))')*epsilon;
    R = R * Exp(delta_R);
    if norm(delta_R)<1e-6
        break
    end
end
R




%% Example 2
clc,clear
% target
f_T = [0,0,1]';

% function
f=@(r)expm(axis2skew (r))*[0,1,0]'; 

% jacobian 
J_Rk=@(r)(eye(3)-(1-cos(norm(r)))/(r'*r)*axis2skew(r)+(norm(r)-sin(norm(r)))/(norm(r)^3)*(axis2skew(r)^2)) ;
J_k = @(r)J_Rk(r)*[0,1,0]';

R = eye(3);

for i =1:100
    tau = skew2axis(logm(R));
    epsilon = f(tau)-f_T;
    delta_R = -((J_k(tau)/J_Rk(tau))'*(J_k(tau)/J_Rk(tau)))\((J_k(tau)/J_Rk(tau))')*epsilon;
    R = R * delta_R;
    if norm(delta_R)<1e-6
        break
    end
end
R

%% Example 3
clc,clear
a=[0 0 1]';
b=[1,0,0]';

R_true=angvec2r(pi/2,[0,1,0]);
R_true*a-b

lb=[0 pi/4 0];
ub=[0 pi*3/4 0];

x=lsqlin(-axis2skew(a),b-a,[],[],[],[],lb,ub)

-axis2skew(a)*x -(b-a)

(eye(3)+axis2skew([0 pi/2 0]))

axis2skew([0 pi/2 0])^2


%% 
clc,clear
syms r1 r2 r3 b1 b2 b3
r=[r1;r2;r3];
b=[b1;b2;b3];
result=axis2skew(r)*b








