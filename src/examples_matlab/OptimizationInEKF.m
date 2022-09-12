%% Load data 
clc,clear
load("data\a_meas.mat");
load("data\w_meas.mat");
load("data\Y_meas.mat");

load("data\X_true.mat");
load("data\bg_true.mat");
load("data\ba_true.mat");
addpath('functionsformytest\')
addpath('MyFunction\')
addpath('Optimization\')

N = length(Y_meas);
Delta_t = 0.001;

%% Cov 
% cov_w_X = blkdiag(0.001*eye(3),0.001*eye(3),zeros(3),0.001*eye(3));
cov_w_X = blkdiag(0.001*eye(3),0.001*eye(3),0.001*eye(3),0.001*eye(3));
% 倒数第二个block 设置为0的话会 MAP 矩阵求导就g了 先这样进行一个比较
cov_bg_all = blkdiag(1e-5*eye(3),1e-5*eye(3));
cov_v_Y = 0.01*eye(3);
M_b = [Delta_t*eye(6);zeros(6)];
H=[zeros(3,6),-eye(3),eye(3)];
d_const = [0 0 0 0 1 -1]';
F = eye(12);
gravity = [ 0 0 -9.81]';
F(4:6,1:3) = axis2skew(Delta_t * gravity);
F(7:9,1:3) = axis2skew(0.5 * Delta_t * Delta_t * gravity);
F(7:9,4:6) = Delta_t * eye(3);

%% InEKF
close all
X_hat_InEKF = zeros(6,6,N);

for i =1:N
    % prediction
    if i==1
        X_k_check = X_true(:,:,1);
        Pk = 0.001*eye(12);
    else
        X_k_check = InEKFprediction(X_hat_InEKF(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1));
        Pk = F*Pk*F'+ cov_w_X;
    end
    % correction
    flag_correction = 1;
    if flag_correction == 1
        H=[zeros(3,6),-eye(3),eye(3)];  % 3*12
        K = Pk*H'/(H*Pk*H'+X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)'); % 12*3
        z = (X_k_check*Y_meas(:,i) - d_const);
        innov = K*z(1:3);
        X_hat_InEKF(:,:,i) = GeneralExp(innov) * X_k_check;
        Pk = (eye(12)-K*H)*Pk;
    else
        X_hat_InEKF(:,:,i) = X_k_check;
    end
end

% position
Plotposition(X_hat_InEKF,X_true,"InEKF")
%% Optimization InEKF
X_hat_MAPInEKF = zeros(6,6,N);

D = [-eye(12),zeros(12);F,-eye(12);zeros(3,12),H];

for i =1:N
    if i ==1 
        X_k_check = X_true(:,:,1);
        D0 = [-eye(12);H];
        W0 = blkdiag(0.001*eye(12),X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
        z = (X_k_check*Y_meas(:,i) - d_const);
        u = [zeros(12,1);z(1:3)];
        sol = (D0'/W0*D0)\D0'/W0*u;
        X_hat_MAPInEKF(:,:,i) = GeneralExp(sol)*X_k_check;
        Pk_inv = (D0'/W0*D0);
    else
        X_k_check = InEKFprediction(X_hat_MAPInEKF(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1));
        W = blkdiag(inv(Pk_inv),cov_w_X,X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
        z = (X_k_check*Y_meas(:,i) - d_const);
        u = [zeros(12,1);zeros(12,1);z(1:3)];
        [sol,Pk_inv,~]=Schur_compute((D'/W*D),D'/W*u);
%         sol = (D'/W*D)\D'/W*u;
        X_hat_MAPInEKF(:,:,i) = GeneralExp(sol)*X_k_check;
    end
end
% position
Plotposition(X_hat_MAPInEKF,X_true,"MAP InEKF")
%% Imperfect InEKF
X_hat_ImInEKF = zeros(6,6,N);
bias_hat = zeros(6,N);

for i =1:N
    % prediction
    if i==1
        X_k_check = X_true(:,:,1);
        bg_bar = zeros(3,1);
        ba_bar = zeros(3,1);
        Pk = 0.001*eye(18);
    else
        X_k_check = InEKFprediction(X_hat_ImInEKF(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1),bg_bar,ba_bar);
        F_ImInEKF = blkdiag(F,eye(6));
        F_ImInEKF(1:12,13:end) = -Adjoint(X_k_check) * M_b;
        Pk = F_ImInEKF*Pk*F_ImInEKF'+ blkdiag(cov_w_X,cov_bg_all);
    end
    % correction
    flag_correction = 1;
    if flag_correction == 1
        H_ImInEKF=[H,zeros(3,6)];  % 3*12
        K = Pk*H_ImInEKF'/(H_ImInEKF*Pk*H_ImInEKF'+X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)'); % 12*3
        z = (X_k_check*Y_meas(:,i) - d_const);
        innov = K*z(1:3);
        X_hat_ImInEKF(:,:,i) = GeneralExp(innov(1:12)) * X_k_check;
        bg_bar = bg_bar + innov(13:15);
        ba_bar = ba_bar + innov(16:end);
        bias_hat(1:3,i) = bg_bar;
        bias_hat(4:6,i) = ba_bar;
        Pk = (eye(18)-K*H_ImInEKF)*Pk;
    else
        X_hat_ImInEKF(:,:,i) = X_k_check;
        bias_hat(1:3,i) = bg_bar;
        bias_hat(4:6,i) = ba_bar;
    end

end
% position
Plotposition(X_hat_ImInEKF,X_true,"Imperfect InEKF")
% bias
Plotbias(bias_hat,bg_true,ba_true)

%% Optimization imperfect InEKF
X_hat_MAPImInEKF = zeros(6,6,N);
bias_hat_MAPImInEKF = zeros(6,N);
X_hat_ImInEKFDebug = zeros(6,6,N);

for i =1:N
    if i ==1 
        X_k_check = X_true(:,:,1);
        bg_check = zeros(3,1);
        ba_check = zeros(3,1);
        D0 = [-eye(18);H,zeros(3,6)];
        W0 = blkdiag(0.001*eye(18),X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
        z = (X_k_check*Y_meas(:,i) - d_const);
        u = [zeros(18,1);z(1:3)];
        sol = (D0'/W0*D0)\D0'/W0*u;
        X_hat_MAPInEKF(:,:,i) = GeneralExp(sol(1:12))*X_k_check;
        bias_hat_MAPImInEKF(1:3,i) = bg_check + sol(13:15);
        bias_hat_MAPImInEKF(4:6,i) = ba_check + sol(16:end);
        Pk_inv = (D0'/W0*D0);
        % debug
        Pk = inv(Pk_inv);
        X_hat_ImInEKFDebug(:,:,i) = X_hat_MAPInEKF(:,:,i);
        bg_bar = sol(13:15);
        ba_bar = sol(16:end);
    else
        X_k_check = InEKFprediction(X_hat_MAPInEKF(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1),...
            bias_hat_MAPImInEKF(1:3,i-1),bias_hat_MAPImInEKF(4:6,i-1));
        bg_check = bias_hat_MAPImInEKF(1:3,i-1);
        ba_check = bias_hat_MAPImInEKF(4:6,i-1);
        W = blkdiag(inv(Pk_inv),cov_w_X,cov_bg_all,X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
        z = (X_k_check*Y_meas(:,i) - d_const);
        u = [zeros(18,1);zeros(18,1);z(1:3)];
        D = -[eye(18),zeros(18);-F,-Adjoint(X_k_check) * M_b,eye(12),zeros(12,6);
            zeros(6,12),-eye(6),zeros(6,12),eye(6);zeros(3,18),-H,zeros(3,6)];
        [sol,Pk_inv,~]=Schur_compute((D'/W*D),D'/W*u);
        X_hat_MAPInEKF(:,:,i) = GeneralExp(sol(1:12))*X_k_check;
        bias_hat_MAPImInEKF(:,i) = [bg_check;ba_check] + sol(13:end);

        % debug
        X_k_check = InEKFprediction(X_hat_ImInEKFDebug(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1),...
            bg_bar,ba_bar);
%         X_k_check = InEKFprediction(X_hat_ImInEKFDebug(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1),...
%             bias_hat_MAPImInEKF(1:3,i-1),bias_hat_MAPImInEKF(4:6,i-1));
        F_ImKF = blkdiag(F,eye(6));
        F_ImKF(1:12,13:end) = -Adjoint(X_k_check) * M_b;
        Pk = F_ImKF*Pk*F_ImKF'+ blkdiag(cov_w_X,cov_bg_all);
        H_IMKF=[zeros(3,6),-eye(3),eye(3),zeros(3,6)];  % 3*12
        K = Pk*H_IMKF'/(H_IMKF*Pk*H_IMKF'+X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)'); % 12*3
        z = (X_k_check*Y_meas(:,i) - d_const);
        innov = K*z(1:3);
        X_hat_ImInEKFDebug(:,:,i) = GeneralExp(innov(1:12)) * X_k_check;
        bg_bar = bg_bar + innov(13:15);
        ba_bar = ba_bar + innov(16:end);
        bias_hat(1:3,i) = bg_bar;
        bias_hat(4:6,i) = ba_bar;
        Pk = (eye(18)-K*H_IMKF)*Pk;
%         bg_bar - bias_hat_MAPImInEKF(1:3,i)
%         ba_bar - bias_hat_MAPImInEKF(4:6,i)
    end
end
% position
Plotposition(X_hat_MAPInEKF,X_true,"MAP Imperfect InEKF")
Plotbias(bias_hat_MAPImInEKF,bg_true,ba_true)
Plotposition(X_hat_ImInEKFDebug,X_true,"MAP Imperfect InEKF Debug")


%% Alternating
% plot 
figure 
p_plot = cell(1,3);
p_plot_true = cell(1,3);
for i =1:3
    subplot(3,1,i);
    p_plot{i} = animatedline;
    p_plot_true{i} = animatedline;
end

%
X_hat_Alternating = zeros(6,6,N);
bias_hat_Alternating = zeros(6,N);

D = [-eye(12),zeros(12);zeros(6,24);F,-eye(12);zeros(3,12),H];

maxiter = 1000;
for i =1:1500
    if i ==1 
        X_k_check = X_true(:,:,1);
        bg_check = zeros(3,1);
        ba_check = zeros(3,1);
        D0 = [-eye(18);H,zeros(3,6)];
        W0 = blkdiag(0.001*eye(18),X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
        z = (X_k_check*Y_meas(:,i) - d_const);
        u = [zeros(18,1);z(1:3)];
        sol = (D0'/W0*D0)\D0'/W0*u;
        X_hat_Alternating(:,:,i) = GeneralExp(sol(1:12))*X_k_check;
        bias_hat_Alternating(1:3,i) = bg_check + sol(13:15);
        bias_hat_Alternating(4:6,i) = ba_check + sol(16:end);
        Pk_inv = (D0'/W0*D0);
        % debug
        Pk = inv(Pk_inv);
    else
        %init
        buffer_last = zeros(24+12,1);
        delta_bg_km1 = zeros(3,1);
        delta_ba_km1 = zeros(3,1);
        bias_k_check = bias_hat_Alternating(:,i-1);
        for k = 1:maxiter
            % fix bias
            X_k_check = InEKFprediction(X_hat_Alternating(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1),...
                bias_hat_Alternating(1:3,i-1) + delta_bg_km1,bias_hat_Alternating(4:6,i-1)+delta_ba_km1);
%             W = blkdiag(inv(Pk_inv),cov_w_X,X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
            W = blkdiag(Pk,cov_w_X,X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
           
            z = (X_k_check*Y_meas(:,i) - d_const);
            u = [zeros(12,1);delta_bg_km1;delta_ba_km1;zeros(12,1);z(1:3)];
%             [~,~,~]=Schur_compute((D'/W*D),D'/W*u);
            sol_all = (D'/W*D)\D'/W*u;
            X_k_next = GeneralExp(sol_all(13:24))*X_k_check;

            % fix states
            z_forbias = [sol_all(1:12);zeros(6,1);sol_all(13:24);zeros(6,1)];
%             W_forbias = blkdiag(inv(Pk_inv),cov_bg_all,cov_w_X);
            W_forbias = blkdiag(Pk,cov_bg_all,cov_w_X);
            D_forbias = [zeros(12,12);-eye(6),zeros(6);Adjoint(X_k_next)*M_b,zeros(12,6);eye(6),-eye(6)];
            sol_bias = (D_forbias'/W_forbias*D_forbias)\D_forbias'/W_forbias*z_forbias;
            delta_bg_km1 =  sol_bias(1:3);
            delta_ba_km1 =  sol_bias(4:6);
            bias_k_hat = bias_k_check + sol_bias(7:end);

            if k>2 && norm([sol_all;sol_bias]-buffer_last)<1e-3
                disp("alternating converge at "+num2str(k)+"-th iteration")
                break;
            end
            if k == maxiter
                disp("maximum iteration! at "+num2str(i))
%                 pause
            end
            buffer_last(1:24) = sol_all;
            buffer_last(25:end) = sol_bias;
        end
        % update estimation 
        X_hat_Alternating(:,:,i) = X_k_next;
        bias_hat_Alternating(:,i) = bias_k_hat;
        %update covirance
        D_cov = [eye(18),zeros(18);
                 -F,-Adjoint(X_k_check) * M_b,eye(12),zeros(12,6);
                 zeros(6,12),-eye(6),zeros(6,12),eye(6);
                 zeros(3,12+6),-H,zeros(3,6)];
        W_cov = blkdiag(inv(Pk_inv),cov_w_X,cov_bg_all,X_k_next(1:3,1:3)*cov_v_Y*X_k_next(1:3,1:3)');
        P_all = D_cov' /W_cov *D_cov; 
        [~,Pk_inv,~]=Schur_compute(P_all,zeros(size(P_all,1),1));

        F_ImKF = blkdiag(F,eye(6));
        F_ImKF(1:12,13:end) = -Adjoint(X_k_check) * M_b;
        Pk = F_ImKF*Pk*F_ImKF'+ blkdiag(cov_w_X,cov_bg_all);
        H_IMKF=[zeros(3,6),-eye(3),eye(3),zeros(3,6)];  % 3*12
        K = Pk*H_IMKF'/(H_IMKF*Pk*H_IMKF'+X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)'); % 12*3
%         z = (X_k_check*Y_meas(:,i) - d_const);
%         innov = K*z(1:3);
%         X_hat_ImInEKFDebug(:,:,i) = GeneralExp(innov(1:12)) * X_k_check;
%         bg_bar = bg_bar + innov(13:15);
%         ba_bar = ba_bar + innov(16:end);
%         bias_hat(1:3,i) = bg_bar;
%         bias_hat(4:6,i) = ba_bar;
        Pk = (eye(18)-K*H_IMKF)*Pk;
%         if abs(1-norm(Pk_inv*Pk_KF))>1e-6
%             error("Pk error")
%         end
    end
    % plot
    for i_plot = 1:3
        addpoints(p_plot{i_plot},i,X_hat_Alternating(i_plot,5,i));
        addpoints(p_plot_true{i_plot},i,X_true(i_plot,5,i));
    end
    drawnow
end


% Plotposition(X_hat_Alternating,X_true,"Alternating")
% Plotbias(bias_hat_Alternating,bg_true,ba_true)


%% My iterative algorithm
X_hat_Myiter = zeros(6,6,N);
bias_hat_Myiter = zeros(6,N);

max_iter = 50;
for i = 1:N
    if i == 1
        %init
        X_0_check = X_true(:,:,1);
        b_0_check = zeros(6,1);
        X_0_bar = X_0_check;
        b_0_bar = b_0_check;
        for k = 1:max_iter
            D = [-eye(18);-H,zeros(3,6)];
            W = blkdiag(0.001*eye(18),X_0_bar(1:3,1:3) * cov_v_Y * X_0_bar(1:3,1:3)');
            z = X_0_bar*Y_meas(:,i)-d_const;
            u = [GeneralLog(X_0_bar*Liegroup_inverse(X_0_check));b_0_bar-b_0_check;z(1:3)];
            sol = (D'/W*D)\D'/W*u;
            %update
            X_0_bar = GeneralExp(sol(1:12)) * X_0_bar;
            b_0_bar = b_0_bar + sol(13:end);
            if norm(sol)<1e-6
                break;
            end
        end
        X_hat_Myiter(:,:,i) = X_0_bar;
        bias_hat_Myiter(:,i) = b_0_bar;
    else
        %init
        

    end
end










