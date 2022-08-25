clc,clear
close all
addpath('MyFunction\')
addpath(genpath('forward_kinematics'))
addpath('Optimization\')

load data\measurements\angular_velocity.mat
load data\measurements\linear_acceleration.mat
load data\measurements\contact.mat
load data\measurements\encoders.mat
load data\measurements\X_init.mat

t = angular_velocity.Time;
angular_mat_ori = angular_velocity.Data;
acc_mat_ori = linear_acceleration.Data;
joint_meas=encoders.Data;
contact_mat = contact.Data;

% ture state
load data\ground_truth\orientation.mat
load data\ground_truth\position.mat
load data\ground_truth\velocity.mat

%% MAIN
N0=10;
N=length(t);    %total length

gravity = -9.81;

w_meas = angular_mat_ori';
a_meas = acc_mat_ori';
encoder_meas = joint_meas';


cov_Ri = 0.001*eye(3);
cov_vi = 0.001*eye(3);
cov_pi = 0.001*eye(3);
cov_bgi = 0.0001*eye(3);
cov_bai = 0.0001*eye(3);
cov_angular_vel = 0.001*eye(3);
cov_acc = 0.001*eye(3);
cov_bgij = 1e-6*eye(3);  
cov_baij = 1e-6*eye(3); 


% Preintegration estimation
Rj_Estimation=zeros(3,3);
vj_Estimation=zeros(3,1);
pj_Estimation=zeros(3,1);
bgj_Estimation=zeros(3,1);
baj_Estimation=zeros(3,1);

Delta_t = 0.0005;

for num_loop = 1:17
    i = (num_loop-1)*1000+N0;
    j = N0+1000*num_loop;

    if num_loop == 1
        check_Ri=orientation.Data(:,:,N0);
        check_vi=velocity.Data(N0,:)';
        check_pi=position.Data(N0,:)';
        check_bgi = zeros(3,1);
        check_bai = zeros(3,1);
        bar_bg = zeros(3,1);
        bar_ba = zeros(3,1);
        
    else
        check_Ri=Rj_Estimation;
        check_vi=vj_Estimation;
        check_pi=pj_Estimation;
        check_bgi = bgj_Estimation;
        check_bai = baj_Estimation;
        bar_bg = check_bgi;
        bar_ba = check_bai;
    end

    % Leg-odometry
    v_tilde_leg = zeros(3,j-i);
    for k = i:j-1
        if contact_mat(k,1)==1 && contact_mat(k,2)== 1
            v_tilde_leg ;
        elseif contact_mat(k,1)==1

        elseif contact_mat(k,2)==1

        end
    end

    % GN initialization
    Ri = check_Ri;
    vi = check_vi;
    pi = check_pi;

    [Bar_Delta_R,Bar_Delta_v,Bar_delta_p,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba]...
        = DeltabarRvpij(i,j,w_meas,bar_bg,a_meas,bar_ba,Delta_t);

    Rj = Ri;
    vj = vi;
    pj = pi;
    for k=i:j-1
        pj = pj + vj*Delta_t + 0.5*(gravity+Rj*(a_meas(:,k)-bar_ba))*Delta_t*Delta_t;
        vj = vj + (gravity+Rj*(a_meas(:,k)-bar_ba))*Delta_t;
        Rj = Rj*Exp((w_meas(:,k)-bar_bg)*Delta_t);
    end
    
    delta_bg = zeros(3,1);
    delta_ba = zeros(3,1);
    maximal_steps = 10;

    for iter_GN = 1:maximal_steps
        %----------------- residual funciton ------------------
        % Prior residual
        r_Ri = Log(check_Ri'*Ri);
        r_vi = vi - check_vi;
        r_pi = pi - check_pi;
        r_bgi = delta_bg;
        r_bai = delta_ba;
        % IMU residual
        r_deltaRij = Log((Bar_Delta_R*Exp(par_R_par_bg*delta_bg))'*Ri'*Rj);
        r_deltavij = Ri'*(vj-vi-gravity*(j-i)*Delta_t)-(Bar_Delta_v+par_v_par_bg*delta_bg+par_v_par_ba*delta_ba);
        r_deltapij = Ri'*(pj-pi-vi*(j-i)*Delta_t-0.5*(j-i)*(j-i)*Delta_t*Delta_t)-...
            (Bar_delta_p+par_p_par_bg*delta_bg+par_p_par_ba*delta_ba);
        % r_bgij = bgj-bgi;  这里先不管 bgj 和 baj
        % r_baij = baj-bai;
        % error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij;r_bgij;r_baij];

        % Leg residual
        % r_dektalij = Ri *(pj-pi) -  Deltatilde_l_ij(i,j,w_meas,bar_bg,delta_bg,Delta_t,tilde_v_k);

        error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij];
        %------------------------ end -------------------------

        %--------------------- covariance ---------------------
        cov_IMU=computeIMUcovariance(i,j,w_meas,bar_bg,delta_bg,a_meas,bar_ba,delta_ba,Delta_t,cov_angular_vel,cov_acc);
        % COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU,cov_bgij,cov_baij);
        COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU);
        %------------------------ end -------------------------

        %------------Jacobian of residual funciton-------------
        % order of jacobian: phii vi pi bgi bai phij vj pj
        par_Ri_par_x = [J_R_inv(r_Ri),zeros(3,21)];
        par_vi_par_x = [zeros(3),eye(3),zeros(3,18)];
        par_pi_par_x = [zeros(3,6),Ri,zeros(3,15)];
        par_bgi_par_x = [zeros(3,9),eye(3),zeros(3,12)];
        par_bai_par_x = [zeros(3,12),eye(3),zeros(3,9)];
        [JacobianR,Jacobianv,Jacobianp] = JacobianofIMUresidualfunc(i,j,Ri,Rj,vi,vj,pi,pj,r_deltaRij,delta_bg,Delta_t,...
            par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba,gravity);
        H = -[par_Ri_par_x;par_vi_par_x;par_pi_par_x;par_bgi_par_x;par_bai_par_x;JacobianR;Jacobianv;Jacobianp];
        %------------------------ end -------------------------

        %---------------- compute the optimal -----------------
        % order: phii vi pi bgi bai phij vj pj
        optimal_sol = (H'/COV_ALL*H)\H'/COV_ALL*error_total;
        %------------------------ end -------------------------

        %------------------------ update ----------------------
        Ri = Ri * Exp(optimal_sol(1:3));
        vi = vi + optimal_sol(4:6);
        pi = pi +Ri * optimal_sol(7:9);
        delta_bg = delta_bg + optimal_sol(10:12);
        delta_ba = delta_ba + optimal_sol(13:15);
        Rj = Rj * Exp(optimal_sol(16:18));
        vj = vj + optimal_sol(19:21);
        pj = pj + Rj * optimal_sol(22:24);
        %------------------------ end -------------------------

        if norm(optimal_sol)<1e-6
            disp("Gaussian-Newton algorithm stoped at "+num2str(iter_GN)+"-th iteration")
            break
        end
    end
  
    % Next loop
    Rj_Estimation = Rj;
    vj_Estimation = vj;
    pj_Estimation = pj;
    bgj_Estimation = bar_bg + delta_bg;
    baj_Estimation = bar_ba + delta_ba;
    
end

%% Single loop test
gravity = -9.81;

w_meas = angular_mat_ori(10:1000,:)';
a_meas = acc_mat_ori(10:1000,:)';
encoder_meas = joint_meas(10:1000,:)';

% initialization
i=1;
j=902;
bar_bg = zeros(3,1);
bar_ba = zeros(3,1);
Delta_t = 0.0005;

check_Ri=eye(3);
check_vi=zeros(3,1);
check_pi=zeros(3,1);
check_bgi = zeros(3,1);
check_bai = zeros(3,1);

cov_Ri = 0.001*eye(3);
cov_vi = 0.001*eye(3);
cov_pi = 0.001*eye(3);
cov_bgi = 0.0001*eye(3);
cov_bai = 0.0001*eye(3);
cov_angular_vel = 0.001*eye(3);
cov_acc = 0.001*eye(3);
cov_bgij = 1e-6*eye(3);  
cov_baij = 1e-6*eye(3); 

Ri = Exp([0.1,0.2,0.1]');
vi = 0.01*ones(3,1);
pi = 0.01*ones(3,1);
% bgi = zeros(3,1);
% bai = zeros(3,1);

[Bar_Delta_R,Bar_Delta_v,Bar_delta_p,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba]...
    = DeltabarRvpij(i,j,w_meas,bar_bg,a_meas,bar_ba,Delta_t);
% Rj = Ri * Bar_Delta_R;
Rj = Ri;
vj = vi;
pj = pi;
for k=i:j-1
    pj = pj + vj*Delta_t + 0.5*(gravity+Rj*(a_meas(:,k)-bar_ba))*Delta_t*Delta_t;
    vj = vj + (gravity+Rj*(a_meas(:,k)-bar_ba))*Delta_t;
    Rj = Rj*Exp((w_meas(:,k)-bar_bg)*Delta_t);
end
% bgj = bgi;
% baj = bai;

delta_bg = zeros(3,1);
delta_ba = zeros(3,1);
maximal_steps = 10;

for iter_GN = 1:maximal_steps
%----------------- residual funciton ------------------
% Prior residual
r_Ri = Log(check_Ri'*Ri);
r_vi = vi - check_vi;
r_pi = pi - check_pi;
r_bgi = delta_bg;
r_bai = delta_ba;

% IMU residual
r_deltaRij = Log((Bar_Delta_R*Exp(par_R_par_bg*delta_bg))'*Ri'*Rj);
r_deltavij = Ri'*(vj-vi-gravity*(j-i)*Delta_t)-(Bar_Delta_v+par_v_par_bg*delta_bg+par_v_par_ba*delta_ba);
r_deltapij = Ri'*(pj-pi-vi*(j-i)*Delta_t-0.5*(j-i)*(j-i)*Delta_t*Delta_t)-...
    (Bar_delta_p+par_p_par_bg*delta_bg+par_p_par_ba*delta_ba);
% r_bgij = bgj-bgi;  这里先不管 bgj 和 baj
% r_baij = baj-bai;
% error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij;r_bgij;r_baij];

% Leg residual
% r_dektalij = Ri *(pj-pi) -  Deltatilde_l_ij(i,j,w_meas,bar_bg,delta_bg,Delta_t,tilde_v_k);

error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij];
%------------------------ end -------------------------

%--------------------- covariance ---------------------
cov_IMU=computeIMUcovariance(i,j,w_meas,bar_bg,delta_bg,a_meas,bar_ba,delta_ba,Delta_t,cov_angular_vel,cov_acc);
% COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU,cov_bgij,cov_baij);
COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU);
%------------------------ end -------------------------

%------------Jacobian of residual funciton-------------
% order of jacobian: phii vi pi bgi bai phij vj pj
par_Ri_par_x = [J_R_inv(r_Ri),zeros(3,21)];
par_vi_par_x = [zeros(3),eye(3),zeros(3,18)];
par_pi_par_x = [zeros(3,6),Ri,zeros(3,15)];
par_bgi_par_x = [zeros(3,9),eye(3),zeros(3,12)];
par_bai_par_x = [zeros(3,12),eye(3),zeros(3,9)];
[JacobianR,Jacobianv,Jacobianp] = JacobianofIMUresidualfunc(i,j,Ri,Rj,vi,vj,pi,pj,r_deltaRij,delta_bg,Delta_t,...
    par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba,gravity);
H = -[par_Ri_par_x;par_vi_par_x;par_pi_par_x;par_bgi_par_x;par_bai_par_x;JacobianR;Jacobianv;Jacobianp];
%------------------------ end -------------------------

%---------------- compute the optimal -----------------
% order: phii vi pi bgi bai phij vj pj
optimal_sol = (H'/COV_ALL*H)\H'/COV_ALL*error_total;
%------------------------ end -------------------------

%------------------------ update ----------------------
Ri = Ri * Exp(optimal_sol(1:3));
vi = vi + optimal_sol(4:6);
pi = pi +Ri * optimal_sol(7:9);
delta_bg = delta_bg + optimal_sol(10:12);
delta_ba = delta_ba + optimal_sol(13:15);
Rj = Rj * Exp(optimal_sol(16:18));
vj = vj + optimal_sol(19:21);
pj = pj + Rj * optimal_sol(22:24);
%------------------------ end -------------------------

    if norm(optimal_sol)<1e-6
        disp("Gaussian-Newton algorithm stoped at "+num2str(iter_GN)+"-th iteration")
        break
    end
end
Ri
Rj
% optimization results

%% Unit Test
w_meas = angular_mat_ori';
a_meas = acc_mat_ori';
bar_bg = zeros(3,1);
bar_ba = zeros(3,1);
Delta_t = 0.0005;

delta_bg = 0.01*ones(3,1);
delta_ba = 0.01*ones(3,1);


[Bar_Delta_R,Bar_Delta_v,Bar_delta_p,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba]...
    = DeltabarRvpij(100,1000,w_meas,bar_bg,a_meas,bar_ba,Delta_t);

%no approximation
w_meas = angular_mat_ori'-delta_bg;
a_meas = acc_mat_ori'-delta_ba;
[Dleta_tilde_R_true,Dleta_tilde_v_true,Dleta_tilde_p_true,~,~,~,~,~]...
    = DeltabarRvpij(100,1000,w_meas,bar_bg,a_meas,bar_ba,Delta_t);


Bar_Delta_R;
norm(Bar_Delta_R-Dleta_tilde_R_true,'fro')
Dleta_tilde_R = Bar_Delta_R*Exp(par_R_par_bg * delta_bg);
norm(Dleta_tilde_R-Dleta_tilde_R_true,'fro')
Dleta_tilde_R_true;

Bar_Delta_v;
norm(Bar_Delta_v-Dleta_tilde_v_true)
Dleta_tilde_v = Bar_Delta_v + par_v_par_bg *delta_bg + par_v_par_ba *delta_ba;
norm(Dleta_tilde_v-Dleta_tilde_v_true)
Dleta_tilde_v_true;

Bar_delta_p;
norm(Bar_delta_p-Dleta_tilde_p_true)
Dleta_tilde_p = Bar_delta_p + par_p_par_bg *delta_bg + par_p_par_ba *delta_ba;
norm(Dleta_tilde_p-Dleta_tilde_p_true)
Dleta_tilde_p_true;


%% TEST

w_meas = angular_mat_ori(100:1000,:)';
a_meas = acc_mat_ori(100:1000,:)';
R0 = eye(3);
v0 = zeros(3,1);
p0 = zeros(3,1);
deltat=0.0005;
bg = 0.001*ones(3,1);
ba = 0.001*ones(3,1);

[R,v,p]=generateIMUdata(w_meas,bg,a_meas,ba,R0,v0,p0,deltat,gravity);
tic
[Dleta_tilde_R,Dleta_tilde_v,Dleta_tilde_p,~,~,~,~,~]...
    = DeltabarRvpij(1,902,w_meas,bg,a_meas,ba,deltat);
disp("true Delta R:")
R(:,:,1)'*R(:,:,end)
Dleta_tilde_R

toc









