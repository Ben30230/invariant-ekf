%% optimization on SO(3) & KF test
addpath('functionsformytest\')
addpath('MyFunction\')
addpath('Optimization\')

%% For optimization
clc,clear
a=[0 0 1]';
b=[1 0 0]';
R_true = angvec2r(pi/2,[0 1 0])
% Rop = Exp(rand(3,1));
Rop = eye(3);
alpha = 0.1;
for i=1:1000
    u = (Rop*a-b);
    J = -axis2skew(Rop*a);
%     rank(J'*J)
%     phi = J'*J\J'*u;
    phi = -alpha *(-a'*Rop'*axis2skew(b))';
    Rop = Exp(phi)*Rop;
    if norm(phi)<1e-6
        break;
    end
end
Rop

%%
clc,clear
a=[0 0 1]';
b=[1 0 0]';
Rop = eye(3);
alpha = 0.1;
for i=1:1000
    delta = -(Rop*a-b)'*axis2skew(Rop*a);
    if abs(delta-(-a'*Rop'*axis2skew(b))')>1e-6
        error("!")
    end
    phi = -alpha *delta';
    Rop = Exp(phi)*Rop;
    if norm(phi)<1e-6
        break;
    end
end
Rop

%%
clc,clear
a=[0 0 1]';
b=[1 0 0]';
Rop = eye(3);
alpha = 0.1;
for i=1:1000
    u = (Rop*a-b);
    J = -axis2skew(Rop*a);
    phi = lsqr(J'*J,-J'*u);
    Rop = Exp(phi)*Rop;
    if norm(phi)<1e-6
        break;
    end
end
Rop
%% linear KF
clc,clear
hat_x_km1 = [1,2]';
hat_P_km1 = [0.01,0.02;0.02,0.02];
Q = 0.001*eye(2);
R = 0.002*eye(2);
F = rand(2);
C = rand(2);
x_check = F*hat_x_km1;
P_k_check = F*hat_P_km1*F' + Q;
K = P_k_check*C'/(C*P_k_check*C'+R);
x_hat = x_check + K * (3-C*x_check)
Pk = (eye(2)-K*C)*P_k_check

% error dynamics version
z = [zeros(2,1);zeros(2,1);3-C*x_check];
H = [-eye(2),zeros(2);F,-eye(2);zeros(2),C];
W = blkdiag(hat_P_km1,Q,R);
A_mat = H'/W*H;
b_mat = H'/W*z;
Pk_inv=A_mat(3:4,3:4)-A_mat(3:4,1:2)/A_mat(1:2,1:2)*A_mat(1:2,3:4);
u = b_mat(3:4)-A_mat(3:4,1:2)/A_mat(1:2,1:2)*b_mat(1:2);
sol = Pk_inv\u;
% e_star = (H'/W*H)\H'/W*z;
x_check + sol
inv(Pk_inv)

% 结果一样的

%% EKF测试
clc,clear
formatSpec = '% 15.8f';

hat_x_km1 = [1 1 0.1]';
P_km1_hat = [0.01 0.03 0
             0.03 0.01 0
             0    0    0.02];
y=[1.5,3,0.1]';

F_Jacobian = @(x)([1+x(3),0.5,x(1);2,1+x(3),x(2);0,0,1]);
G_Jacobian = @(x)(2*diag([x(1),x(2),x(3)]));

Temp_C = rand(3);
Q = 0.01*(Temp_C*Temp_C');
R = 0.02*eye(3);
check_x_k = f_sys(hat_x_km1);
check_P_k = F_Jacobian(hat_x_km1) * P_km1_hat * F_Jacobian(hat_x_km1)'+Q;


hat_x = check_x_k;
hat_x_last = check_x_k;
for i=1:10
    K = check_P_k*G_Jacobian(hat_x)'/(G_Jacobian(hat_x)*check_P_k*G_Jacobian(hat_x)'+R);
    innovation = K*(y-g_sys(hat_x)-G_Jacobian(hat_x)*(check_x_k-hat_x));
    hat_x = check_x_k + innovation;
    if i==1
        disp("EKF: x_hat = "+num2str(hat_x',formatSpec))
    end
    if norm(hat_x-hat_x_last)<1e-6
        break;
    end
    hat_x_last = hat_x;
end
disp("IEKF: x_hat = "+num2str(hat_x',formatSpec)+"  at "+num2str(i)+"-th iteration")

% MAP
x_km1_op = hat_x_km1;
x_k_op = check_x_k;
W = blkdiag(P_km1_hat,Q,R);
for i = 1:1000   
    H = [eye(3),zeros(3);-F_Jacobian(x_km1_op),eye(3);zeros(3),G_Jacobian(x_k_op)];
    rank(H);
    e = [hat_x_km1-x_km1_op;f_sys(x_km1_op)-x_k_op;y-g_sys(x_k_op)];
    delta_x = (H'/W*H)\H'/W * e;
%     delta_x = - (H'*H)\H' * e;
    x_km1_op = x_km1_op +delta_x(1:3);
    x_k_op=x_k_op+delta_x(4:6);
    if i==1
        disp("MAP at first iteration: x_hat = "+num2str(x_k_op',formatSpec))
    end
    if norm(delta_x)<1e-6
        break;
    end
end
disp("MAP: x_hat = "+num2str(x_k_op',formatSpec)+"  at "+num2str(i)+"-th iteration")

% MAP version 2
x_k_op = check_x_k;
for i=1:50
    H = [eye(3);G_Jacobian(x_k_op)];
    W = blkdiag(check_P_k,R);
    e = [check_x_k-x_k_op;y-g_sys(x_k_op)];
    delta = (H'/W*H)\H'/W*e;
    x_k_op = x_k_op +delta;
    if i==1
        disp("MAP (second version) at first iteration: x_hat = "+num2str(x_k_op',formatSpec))
    end
    if norm(delta)<1e-6
        break;
    end
end
disp("MAP (second version): x_hat = "+num2str(x_k_op',formatSpec)+"  at "+num2str(i)+"-th iteration")







%%
a = pi/2*[0 0 1]';
J_R(a*0.001)





























