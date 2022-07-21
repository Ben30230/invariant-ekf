clc,clear
close all
addpath('MyFunction\')
addpath(genpath('forward_kinematics'))

load data\measurements\angular_velocity.mat
load data\measurements\linear_acceleration.mat
load data\measurements\contact.mat
load data\measurements\encoders.mat
load data\measurements\X_init.mat

t = angular_velocity.Time;
angular_mat = angular_velocity.Data;
acc_mat = linear_acceleration.Data;
joint_meas=encoders.Data;
contact_mat = contact.Data;

% ture state
load data\ground_truth\orientation.mat
load data\ground_truth\position.mat
load data\ground_truth\velocity.mat

%%
N0= 10;
R0=orientation.Data(:,:,N0);
v0=velocity.Data(N0,:)';
p0=position.Data(N0,:)';
dl=R0*p_VectorNav_to_LeftToeBottom(joint_meas(N0,:));
dr=R0*p_VectorNav_to_RightToeBottom(joint_meas(N0,:));
d_all=[dl,dr];
P=0.001*eye(9);
bg = zeros(3,1);
ba = zeros(3,1);
robotstate=RobotState_Bias(R0,v0,p0,d_all,P,0);

N=length(t);
% N=1000;

% estimation
R_Estimation=zeros(3,3,N);
v_Estimation=zeros(3,N);
p_Estimation=zeros(3,N);

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

    robotstate.prediction(angular_mat(i,:),acc_mat(i,:),joint_meas(i,:),0.0005);
    
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
    hold off
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
figure
for i=1:3
    subplot(3,1,i)
    plot(N0:N,p_Estimation(i,N0:N))
    hold on
    plot(N0:N,position.Data(N0:N,i))
    hold off
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

% % SO3
% error_R=zeros(3,N);
% for i=N0:N
%     tem_logcoordinates=logm(orientation.Data(:,:,i)*R_Estimation(:,:,i)');
%     error_R(:,i)=skew2axis(tem_logcoordinates);
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