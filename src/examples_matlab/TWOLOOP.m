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

% encoder_diff
encoder_dot_meas = zeros(14,length(t));
for i = 2: length(t)-1
    encoder_dot_meas(:,i) = (joint_meas(i+1,:)-joint_meas(i-1,:))'/(2*0.0005);
end

%% MAIN
leg_flag = 1 ;

N0=1265;
horizon = 500;
N=length(t);    %total length

gravity = [0,0,-9.81]';

w_meas = angular_mat_ori';
a_meas = acc_mat_ori';
encoder_meas = joint_meas;


cov_Ri = 0.001*eye(3);
cov_vi = 0.001*eye(3);
cov_pi = 0.001*eye(3);
cov_bgi = 0.0001*eye(3);
cov_bai = 0.0001*eye(3);
cov_angular_vel = 0.001*eye(3);
cov_acc = 0.001*eye(3);
cov_legdrift = 0.01*eye(3);
% cov_bgij = 1e-6*eye(3);  
% cov_baij = 1e-6*eye(3); 


% Preintegration estimation
Rj_Estimation=zeros(3,3);
vj_Estimation=zeros(3,1);
pj_Estimation=zeros(3,1);
bgj_Estimation=zeros(3,1);
baj_Estimation=zeros(3,1);

% estimation
R_Estimation=zeros(3,3,N);
v_Estimation=zeros(3,N);
p_Estimation=zeros(3,N);

bg_preinter=zeros(3,N);
ba_preinter=zeros(3,N);

Delta_t = 0.0005;

% Pure IMU
R_only_IMU=zeros(3,3,N);
v_only_IMU=zeros(3,N);
p_only_IMU=zeros(3,N);
R_only_IMU(:,:,N0) = orientation.Data(:,:,N0);
v_only_IMU(:,N0) = velocity.Data(N0,:)';
p_only_IMU(:,N0) = position.Data(N0,:)';
for i = N0:N-1
    R_only_IMU(:,:,i+1) = R_only_IMU(:,:,i)*Exp((w_meas(:,i)-0)*Delta_t);
    v_only_IMU(:,i+1) = v_only_IMU(:,i) + (gravity+R_only_IMU(:,:,i)*(a_meas(:,i)-0))*Delta_t;
    p_only_IMU(:,i+1) = p_only_IMU(:,i) + v_only_IMU(:,i)*Delta_t + 0.5*(gravity+R_only_IMU(:,:,i)*(a_meas(:,i)-0))*Delta_t*Delta_t;
end
% figure
% for i=1:3
%     subplot(3,1,i)
%     plot(N0:N,v_only_IMU(i,N0:N))
% end
% figure
% for i=1:3
%     subplot(3,1,i)
%     plot(N0:N,p_only_IMU(i,N0:N))
% end

%
for num_loop = 1:37
    i = (num_loop-1)*horizon+N0;
    j = N0+horizon*num_loop;

    if num_loop == 1
        check_Ri=orientation.Data(:,:,N0);
        check_vi=velocity.Data(N0,:)';
        check_pi=position.Data(N0,:)';
%         check_bgi = zeros(3,1);
%         check_bai = zeros(3,1);
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
%     v_tilde_leg = zeros(3,j-i);
%     for k = i:j-1
%         if contact_mat(k,1)==1 && contact_mat(k,2)== 1
%             v_tilde_leg ;
%         elseif contact_mat(k,1)==1
% 
%         elseif contact_mat(k,2)==1
% 
%         end
%     end

    % GN initialization
    Ri = check_Ri;
%     Ri =eye(3);
    vi = check_vi;
    pi = check_pi;

%     [Bar_Delta_R,Bar_Delta_v,Bar_delta_p,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba]...
%         = DeltabarRvpij(i,j,w_meas,bar_bg,a_meas,bar_ba,Delta_t);
    [Bar_Delta_R,Bar_Delta_v,Bar_delta_p,Bar_delta_l,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba,par_l_par_bg]...
        = DeltabarRvpij(i,j,w_meas,bar_bg,a_meas,bar_ba,Delta_t,encoder_meas,contact_mat,encoder_dot_meas);

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
    maximal_steps = 50;

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
        r_deltapij = Ri'*(pj-pi-vi*(j-i)*Delta_t-0.5*gravity * (j-i)*(j-i)*Delta_t*Delta_t)-...
            (Bar_delta_p+par_p_par_bg*delta_bg+par_p_par_ba*delta_ba);
        % r_bgij = bgj-bgi;  这里先不管 bgj 和 baj
        % r_baij = baj-bai;
        % error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij;r_bgij;r_baij];

        % Leg residual
        r_dektalij = Ri' *(pj-pi) -  Bar_delta_l + par_l_par_bg * delta_bg; 
%         Bar_delta_l + par_l_par_bg * delta_bg
%         r_dektalij_true = orientation.Data(:,:,i)*(position.Data(j,:)-position.Data(i,:))'
        
%         Ri' *(pj-pi) - (p_VectorNav_to_RightToeBottom(encoder_meas(i,:)) ...
%             - Bar_Delta_R * p_VectorNav_to_RightToeBottom(encoder_meas(j,:)))

        if leg_flag
            error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij;r_dektalij];
        else
            error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij];
        end
        
        %------------------------ end -------------------------

        %--------------------- covariance ---------------------
        [cov_IMU,Sigma_leg_ik]=computeIMUcovariance(i,j,w_meas,bar_bg,delta_bg,a_meas,bar_ba,...
            delta_ba,Delta_t,cov_angular_vel,cov_acc,cov_legdrift,encoder_meas,contact_mat,encoder_dot_meas);
        % COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU,cov_bgij,cov_baij);
        if leg_flag
            COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU,Sigma_leg_ik);
%             COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU,eye(3));
        else
            COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU);
        end
        
        %------------------------ end -------------------------

        %------------Jacobian of residual funciton-------------
        % order of jacobian: phii vi pi bgi bai phij vj pj
        par_Ri_par_x = [J_R_inv(r_Ri),zeros(3,21)];
        par_vi_par_x = [zeros(3),eye(3),zeros(3,18)];
        par_pi_par_x = [zeros(3,6),Ri,zeros(3,15)];
        par_bgi_par_x = [zeros(3,9),eye(3),zeros(3,12)];
        par_bai_par_x = [zeros(3,12),eye(3),zeros(3,9)];
        [JacobianR,Jacobianv,Jacobianp,Jacobianl] = JacobianofIMUresidualfunc(i,j,Ri,Rj,vi,vj,pi,pj,r_deltaRij,delta_bg,Delta_t,...
            par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba,gravity,par_l_par_bg);

        if leg_flag
            H = -[par_Ri_par_x;par_vi_par_x;par_pi_par_x;par_bgi_par_x;par_bai_par_x;JacobianR;Jacobianv;Jacobianp;Jacobianl];
        else
            H = -[par_Ri_par_x;par_vi_par_x;par_pi_par_x;par_bgi_par_x;par_bai_par_x;JacobianR;Jacobianv;Jacobianp];
        end
        
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
        if iter_GN == maximal_steps
            disp("maximal steps!     norm:"+ num2str(norm(optimal_sol)))
        end
    end
  
    % Next loop
    Rj_Estimation = Rj;
    vj_Estimation = vj;
    pj_Estimation = pj;
    bgj_Estimation = bar_bg + delta_bg;
    baj_Estimation = bar_ba + delta_ba;
    
    % record
    R_Estimation(:,:,j) = Rj_Estimation;
    v_Estimation(:,j)=vj_Estimation;
    p_Estimation(:,j)=pj_Estimation;
    bg_preinter(:,i:j-1) = repmat(bgj_Estimation,[1,horizon]) ;
    ba_preinter(:,i:j-1) = repmat(baj_Estimation,[1,horizon]) ;
end

%% PLOT
close all

% velocity
v_preintegration = v_Estimation(:,sum(v_Estimation)~=0);
v_preintegration_index = find(sum(v_Estimation)~=0);
figure
for i=1:3
    subplot(3,1,i)
    plot(v_preintegration_index,v_preintegration(i,:),'o')
    hold on
    plot(N0:N,velocity.Data(N0:N,i))
%     hold on
%     plot(N0:N,v_only_IMU(i,N0:N))
    hold off
%     legend('Estiamtion','GroundTruth','Only IMU')
    legend('Estiamtion','GroundTruth')
    xlabel('time step')
    switch i
        case 1
            ylabel('x(m/s)')
        case 2
            ylabel('y(m/s)')
        case 3
            ylabel('z(m/s)')
    end
end
sgtitle("Estimation of linear velocity")

% position
p_preintegration = p_Estimation(:,sum(p_Estimation)~=0);
p_preintegration_index = find(sum(p_Estimation)~=0);
figure
for i=1:3
    subplot(3,1,i)
    plot(p_preintegration_index,p_preintegration(i,:),'o')
    hold on
    plot(N0:N,position.Data(N0:N,i))
%     hold on
%     plot(N0:N,p_only_IMU(i,N0:N))
    hold off
%     legend('Estiamtion','GroundTruth','Only IMU')
    legend('Estiamtion','GroundTruth')
    xlabel('time step')
    switch i
        case 1
            ylabel('x(m)')
        case 2
            ylabel('y(m)')
        case 3
            ylabel('z(m)')
    end
end
sgtitle("Estimation of linear position")

% % bias
% figure
% for i=1:2
%     for j=1:3
%         subplot(3,2,i+2*(j-1))     
%         if i==1
%             plot(N0:N,bg_Estimation(j,N0:N))
%             title('bg')
%         else
%             plot(N0:N,ba_Estimation(j,N0:N))
%             title('ba')
%         end
%     end
% end
% sgtitle("Estimation of Bias")

% % SO3
% error_R=zeros(3,N);
% for i=N0:N
%     error_R(:,i)=Log(orientation.Data(:,:,i)*R_Estimation(:,:,i)');
% end
% figure
% for i=1:3
%     subplot(3,1,i)
%     plot(N0:N,error_R(i,N0:N))
%     xlabel('time step')
%     switch i
%         case 1
%             ylabel('log-coor(1)')
%         case 2
%             ylabel('log-coor(2)')
%         case 3
%             ylabel('log-coor(3)')
%     end
% end
% sgtitle("Estimation of orientation (log(R_tR'))")

%% InEKF no bias
N0= 1263;
R0=orientation.Data(:,:,N0);
v0=velocity.Data(N0,:)';
p0=position.Data(N0,:)';
dl=R0*p_VectorNav_to_LeftToeBottom(joint_meas(N0,:));
dr=R0*p_VectorNav_to_RightToeBottom(joint_meas(N0,:));
d_all=[dl,dr];
flag_bias = 0;
if flag_bias
    P=0.001*eye(15);
else
    P=0.001*eye(9);
end
bg = zeros(3,1);
ba = zeros(3,1);
robotstate=RobotState_Bias(R0,v0,p0,d_all,P,flag_bias);

N=length(t);
% N=1000;

% bias constant
bg_offset = zeros(3,1);
ba_offset = zeros(3,1);
bg_array = zeros(3,1000);
ba_array = zeros(3,1000);
bg_array_real = zeros(3,1000);
k_b=1;

% estimation
R_Estimation=zeros(3,3,N);
v_Estimation=zeros(3,N);
p_Estimation=zeros(3,N);

bg_Estimation=zeros(3,N);
ba_Estimation=zeros(3,N);

angular_mat = angular_mat_ori;
acc_mat = acc_mat_ori;
for i=N0:N
    if contact_mat(i,1) == 1
        robotstate.legcontact_flag_member(1) = 1;
    else
        robotstate.legcontact_flag_member(1) = 0;
    end
    if contact_mat(i,2) == 1
        robotstate.legcontact_flag_member(2) = 1;
    else
        robotstate.legcontact_flag_member(2) = 0;
    end

    if sum(robotstate.legcontact_flag_member == robotstate.legcontact_last_flag_member) ~=2
        robotstate.ContactUpdata(joint_meas(i,:));
%         disp("switch at i="+num2str(i))
    end


    if i > N0
        robotstate.prediction(angular_mat(i,:),acc_mat(i,:),joint_meas(i,:),0.0005);
    end
    
    %debug 
%     if i==443
%         robotstate.P_member
%         robotstate.X_member
%     end

    if sum(robotstate.legcontact_flag_member) ~= 0
        robotstate.Correction(joint_meas(i,:));
    end

    % record estimation
    R_Estimation(:,:,i) = robotstate.R_member;
    v_Estimation(:,i) = robotstate.v_member;
    p_Estimation(:,i) = robotstate.p_member;
    bg_Estimation(:,i) = robotstate.bg_member;
    ba_Estimation(:,i) = robotstate.ba_member;
end
%% with bias from MAP
N0= 1263;
R0=orientation.Data(:,:,N0);
v0=velocity.Data(N0,:)';
p0=position.Data(N0,:)';
dl=R0*p_VectorNav_to_LeftToeBottom(joint_meas(N0,:));
dr=R0*p_VectorNav_to_RightToeBottom(joint_meas(N0,:));
d_all=[dl,dr];
flag_bias = 0;
if flag_bias
    P=0.001*eye(15);
else
    P=0.001*eye(9);
end
bg = zeros(3,1);
ba = zeros(3,1);
robotstate=RobotState_Bias(R0,v0,p0,d_all,P,flag_bias);

N=length(t);
% N=1000;


% estimation
R_twoloop=zeros(3,3,N);
v_twoloop=zeros(3,N);
p_twoloop=zeros(3,N);

bg_Estimation=zeros(3,N);
ba_Estimation=zeros(3,N);

angular_mat = angular_mat_ori-bg_preinter';
acc_mat = acc_mat_ori-ba_preinter';
for i=N0:N
    if contact_mat(i,1) == 1
        robotstate.legcontact_flag_member(1) = 1;
    else
        robotstate.legcontact_flag_member(1) = 0;
    end
    if contact_mat(i,2) == 1
        robotstate.legcontact_flag_member(2) = 1;
    else
        robotstate.legcontact_flag_member(2) = 0;
    end

    if sum(robotstate.legcontact_flag_member == robotstate.legcontact_last_flag_member) ~=2
        robotstate.ContactUpdata(joint_meas(i,:));
%         disp("switch at i="+num2str(i))
    end


    if i > N0
        robotstate.prediction(angular_mat(i,:),acc_mat(i,:),joint_meas(i,:),0.0005);
    end
    
    %debug 
%     if i==443
%         robotstate.P_member
%         robotstate.X_member
%     end

    if sum(robotstate.legcontact_flag_member) ~= 0
        robotstate.Correction(joint_meas(i,:));
    end

    % record estimation
    R_twoloop(:,:,i) = robotstate.R_member;
    v_twoloop(:,i) = robotstate.v_member;
    p_twoloop(:,i) = robotstate.p_member;
    bg_Estimation(:,i) = robotstate.bg_member;
    ba_Estimation(:,i) = robotstate.ba_member;
end




%% PLOT
close all

% velocity
figure
for i=1:3
    subplot(3,1,i)
    plot(N0:N,v_Estimation(i,N0:N))
    hold on
    plot(N0:N,velocity.Data(N0:N,i))
    hold on
    plot(N0:N,v_twoloop(i,N0:N))
    hold off
    legend('Estiamtion','GroundTruth','Two loop')
    xlabel('time step')
    switch i
        case 1
            ylabel('x(m/s)')
        case 2
            ylabel('y(m/s)')
        case 3
            ylabel('z(m/s)')
    end
end
sgtitle("Estimation of linear velocity")
ATE_velocity_nobias = norm(v_Estimation(:,N0:N)-velocity.Data(N0:N,:)','fro');
ATE_velocity_twoloop = norm(v_twoloop(:,N0:N)-velocity.Data(N0:N,:)','fro');
disp("ATE of position for only InEKF: "+num2str(ATE_velocity_nobias));
disp("ATE of position for two loops: "+num2str(ATE_velocity_twoloop));

% position
figure
for i=1:3
    subplot(3,1,i)
    plot(N0:N,p_Estimation(i,N0:N))
    hold on
    plot(N0:N,position.Data(N0:N,i))
    hold on
    plot(N0:N,p_twoloop(i,N0:N))
    hold off
    legend('Estiamtion','GroundTruth','Two loop')
    xlabel('time step')
    switch i
        case 1
            ylabel('x(m)')
        case 2
            ylabel('y(m)')
        case 3
            ylabel('z(m)')
    end
end
sgtitle("Estimation of linear position")
ATE_position_nobias = norm(p_Estimation(:,N0:N)-position.Data(N0:N,:)','fro');
ATE_position_twoloop = norm(p_twoloop(:,N0:N)-position.Data(N0:N,:)','fro');
disp("ATE of position for only InEKF: "+num2str(ATE_position_nobias));
disp("ATE of position for two loops: "+num2str(ATE_position_twoloop));

% bias
figure
for i=1:2
    for j=1:3
        subplot(3,2,i+2*(j-1))     
        if i==1
            plot(N0:N,bg_Estimation(j,N0:N))
            title('bg')
        else
            plot(N0:N,ba_Estimation(j,N0:N))
            title('ba')
        end
    end
end
sgtitle("Estimation of Bias")

% SO3
error_R=zeros(3,N);
for i=N0:N
    error_R(:,i)=Log(orientation.Data(:,:,i)*R_Estimation(:,:,i)');
end
figure
for i=1:3
    subplot(3,1,i)
    plot(N0:N,error_R(i,N0:N))
    xlabel('time step')
    switch i
        case 1
            ylabel('log-coor(1)')
        case 2
            ylabel('log-coor(2)')
        case 3
            ylabel('log-coor(3)')
    end
end
sgtitle("Estimation of orientation (log(R_tR'))")





% %% Single loop test
% gravity = -9.81;
% 
% w_meas = angular_mat_ori(10:1000,:)';
% a_meas = acc_mat_ori(10:1000,:)';
% encoder_meas = joint_meas(10:1000,:)';
% 
% % initialization
% i=1;
% j=902;
% bar_bg = zeros(3,1);
% bar_ba = zeros(3,1);
% Delta_t = 0.0005;
% 
% check_Ri=eye(3);
% check_vi=zeros(3,1);
% check_pi=zeros(3,1);
% check_bgi = zeros(3,1);
% check_bai = zeros(3,1);
% 
% cov_Ri = 0.001*eye(3);
% cov_vi = 0.001*eye(3);
% cov_pi = 0.001*eye(3);
% cov_bgi = 0.0001*eye(3);
% cov_bai = 0.0001*eye(3);
% cov_angular_vel = 0.001*eye(3);
% cov_acc = 0.001*eye(3);
% cov_bgij = 1e-6*eye(3);  
% cov_baij = 1e-6*eye(3); 
% 
% Ri = Exp([0.1,0.2,0.1]');
% vi = 0.01*ones(3,1);
% pi = 0.01*ones(3,1);
% % bgi = zeros(3,1);
% % bai = zeros(3,1);
% 
% [Bar_Delta_R,Bar_Delta_v,Bar_delta_p,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba]...
%     = DeltabarRvpij(i,j,w_meas,bar_bg,a_meas,bar_ba,Delta_t);
% % Rj = Ri * Bar_Delta_R;
% Rj = Ri;
% vj = vi;
% pj = pi;
% for k=i:j-1
%     pj = pj + vj*Delta_t + 0.5*(gravity+Rj*(a_meas(:,k)-bar_ba))*Delta_t*Delta_t;
%     vj = vj + (gravity+Rj*(a_meas(:,k)-bar_ba))*Delta_t;
%     Rj = Rj*Exp((w_meas(:,k)-bar_bg)*Delta_t);
% end
% % bgj = bgi;
% % baj = bai;
% 
% delta_bg = zeros(3,1);
% delta_ba = zeros(3,1);
% maximal_steps = 10;
% 
% for iter_GN = 1:maximal_steps
% %----------------- residual funciton ------------------
% % Prior residual
% r_Ri = Log(check_Ri'*Ri);
% r_vi = vi - check_vi;
% r_pi = pi - check_pi;
% r_bgi = delta_bg;
% r_bai = delta_ba;
% 
% % IMU residual
% r_deltaRij = Log((Bar_Delta_R*Exp(par_R_par_bg*delta_bg))'*Ri'*Rj);
% r_deltavij = Ri'*(vj-vi-gravity*(j-i)*Delta_t)-(Bar_Delta_v+par_v_par_bg*delta_bg+par_v_par_ba*delta_ba);
% r_deltapij = Ri'*(pj-pi-vi*(j-i)*Delta_t-0.5*(j-i)*(j-i)*Delta_t*Delta_t)-...
%     (Bar_delta_p+par_p_par_bg*delta_bg+par_p_par_ba*delta_ba);
% % r_bgij = bgj-bgi;  这里先不管 bgj 和 baj
% % r_baij = baj-bai;
% % error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij;r_bgij;r_baij];
% 
% % Leg residual
% % r_dektalij = Ri *(pj-pi) -  Deltatilde_l_ij(i,j,w_meas,bar_bg,delta_bg,Delta_t,tilde_v_k);
% 
% error_total = [r_Ri;r_vi;r_pi;r_bgi;r_bai;r_deltaRij;r_deltavij;r_deltapij];
% %------------------------ end -------------------------
% 
% %--------------------- covariance ---------------------
% cov_IMU=computeIMUcovariance(i,j,w_meas,bar_bg,delta_bg,a_meas,bar_ba,delta_ba,Delta_t,cov_angular_vel,cov_acc);
% % COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU,cov_bgij,cov_baij);
% COV_ALL = blkdiag(cov_Ri,cov_vi,cov_pi,cov_bgi,cov_bai,cov_IMU);
% %------------------------ end -------------------------
% 
% %------------Jacobian of residual funciton-------------
% % order of jacobian: phii vi pi bgi bai phij vj pj
% par_Ri_par_x = [J_R_inv(r_Ri),zeros(3,21)];
% par_vi_par_x = [zeros(3),eye(3),zeros(3,18)];
% par_pi_par_x = [zeros(3,6),Ri,zeros(3,15)];
% par_bgi_par_x = [zeros(3,9),eye(3),zeros(3,12)];
% par_bai_par_x = [zeros(3,12),eye(3),zeros(3,9)];
% [JacobianR,Jacobianv,Jacobianp] = JacobianofIMUresidualfunc(i,j,Ri,Rj,vi,vj,pi,pj,r_deltaRij,delta_bg,Delta_t,...
%     par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba,gravity);
% H = -[par_Ri_par_x;par_vi_par_x;par_pi_par_x;par_bgi_par_x;par_bai_par_x;JacobianR;Jacobianv;Jacobianp];
% %------------------------ end -------------------------
% 
% %---------------- compute the optimal -----------------
% % order: phii vi pi bgi bai phij vj pj
% optimal_sol = (H'/COV_ALL*H)\H'/COV_ALL*error_total;
% %------------------------ end -------------------------
% 
% %------------------------ update ----------------------
% Ri = Ri * Exp(optimal_sol(1:3));
% vi = vi + optimal_sol(4:6);
% pi = pi +Ri * optimal_sol(7:9);
% delta_bg = delta_bg + optimal_sol(10:12);
% delta_ba = delta_ba + optimal_sol(13:15);
% Rj = Rj * Exp(optimal_sol(16:18));
% vj = vj + optimal_sol(19:21);
% pj = pj + Rj * optimal_sol(22:24);
% %------------------------ end -------------------------
% 
%     if norm(optimal_sol)<1e-6
%         disp("Gaussian-Newton algorithm stoped at "+num2str(iter_GN)+"-th iteration")
%         break
%     end
% end
% Ri
% Rj
% % optimization results
% 
% %% Unit Test
% w_meas = angular_mat_ori';
% a_meas = acc_mat_ori';
% bar_bg = zeros(3,1);
% bar_ba = zeros(3,1);
% Delta_t = 0.0005;
% 
% delta_bg = 0.01*ones(3,1);
% delta_ba = 0.01*ones(3,1);
% 
% 
% [Bar_Delta_R,Bar_Delta_v,Bar_delta_p,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba]...
%     = DeltabarRvpij(100,1000,w_meas,bar_bg,a_meas,bar_ba,Delta_t);
% 
% %no approximation
% w_meas = angular_mat_ori'-delta_bg;
% a_meas = acc_mat_ori'-delta_ba;
% [Dleta_tilde_R_true,Dleta_tilde_v_true,Dleta_tilde_p_true,~,~,~,~,~]...
%     = DeltabarRvpij(100,1000,w_meas,bar_bg,a_meas,bar_ba,Delta_t);
% 
% 
% Bar_Delta_R;
% norm(Bar_Delta_R-Dleta_tilde_R_true,'fro')
% Dleta_tilde_R = Bar_Delta_R*Exp(par_R_par_bg * delta_bg);
% norm(Dleta_tilde_R-Dleta_tilde_R_true,'fro')
% Dleta_tilde_R_true;
% 
% Bar_Delta_v;
% norm(Bar_Delta_v-Dleta_tilde_v_true)
% Dleta_tilde_v = Bar_Delta_v + par_v_par_bg *delta_bg + par_v_par_ba *delta_ba;
% norm(Dleta_tilde_v-Dleta_tilde_v_true)
% Dleta_tilde_v_true;
% 
% Bar_delta_p;
% norm(Bar_delta_p-Dleta_tilde_p_true)
% Dleta_tilde_p = Bar_delta_p + par_p_par_bg *delta_bg + par_p_par_ba *delta_ba;
% norm(Dleta_tilde_p-Dleta_tilde_p_true)
% Dleta_tilde_p_true;
% 
% 
% %% TEST
% 
% w_meas = angular_mat_ori(100:1000,:)';
% a_meas = acc_mat_ori(100:1000,:)';
% R0 = eye(3);
% v0 = zeros(3,1);
% p0 = zeros(3,1);
% deltat=0.0005;
% bg = 0.001*ones(3,1);
% ba = 0.001*ones(3,1);
% 
% [R,v,p]=generateIMUdata(w_meas,bg,a_meas,ba,R0,v0,p0,deltat,gravity);
% tic
% [Dleta_tilde_R,Dleta_tilde_v,Dleta_tilde_p,~,~,~,~,~]...
%     = DeltabarRvpij(1,902,w_meas,bg,a_meas,ba,deltat);
% disp("true Delta R:")
% R(:,:,1)'*R(:,:,end)
% Dleta_tilde_R
% 
% toc









