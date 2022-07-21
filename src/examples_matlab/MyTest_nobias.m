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


%% 从第N0个数据开始
% X0=X_init.Data(:,:,2);
N0=10;
R0=orientation.Data(:,:,N0);
v0=velocity.Data(N0,:)';
p0=position.Data(N0,:)';
dl=R0*p_VectorNav_to_LeftToeBottom(joint_meas(N0,:));
dr=R0*p_VectorNav_to_RightToeBottom(joint_meas(N0,:));
P=0.001*eye(15);
% P(7:9,7:9)=0;
robotstate=RobotState(R0,v0,p0,dl,dr,P);

N=length(t);
% N=1000;

% estimation
R_Estimation=zeros(3,3,N);
v_Estimation=zeros(3,N);
p_Estimation=zeros(3,N);

for i=N0:N
%     if i==223
%         pause
%     end
    % initialization: leg contacting
    switch robotstate.contact_flag_member
        case 0
            if contact_mat(i,1)==1 && contact_mat(i,2)==1% both
                robotstate.Contactnone2two(joint_meas(i,:));
            elseif contact_mat(i,1)==1 %left only
                robotstate.Contactnone2left(joint_meas(i,:));
            elseif contact_mat(i,2)==1 %right only
                robotstate.Contactnone2right(joint_meas(i,:));
            end
        case 1 % left 
            if contact_mat(i,1)==1 && contact_mat(i,2)==1% both
                robotstate.Contactleft2two(joint_meas(i,:));
            elseif contact_mat(i,2)==1 %right only
                robotstate.Contactleft2right(joint_meas(i,:));
            elseif contact_mat(i,1)~= 1 && contact_mat(i,2)~= 1 %none
                robotstate.Contactleft2none();
            end
        case 2 %right
            if contact_mat(i,1)==1 && contact_mat(i,2)==1% both
                robotstate.Contactright2two(joint_meas(i,:));
            elseif contact_mat(i,1)==1 %left only
                robotstate.Contactright2left(joint_meas(i,:));
            elseif contact_mat(i,1)~= 1 && contact_mat(i,2)~= 1 %none
                robotstate.Contactright2none();
            end
        case 3 % both
            if contact_mat(i,1)~= 1 && contact_mat(i,2)~= 1 %none
                robotstate.Contacttwo2none();
            elseif contact_mat(i,1)== 1 && contact_mat(i,2)~= 1 % left only
                robotstate.Contacttwo2left();
            elseif contact_mat(i,1) ~= 1 && contact_mat(i,2) == 1 % right only
                robotstate.Contacttwo2right();
            end
    end

    % prediction
    robotstate.prediction(angular_mat(i,:),acc_mat(i,:),joint_meas(i,:),0.0005);

    %debug 
%     if i==443
%         robotstate.P_member
%         robotstate.X_member
%     end

    % correction
    if robotstate.contact_flag_member~=0 && i>20
        robotstate.Correction(joint_meas(i,:));
%         disp("pudate!")
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
% for i=2:N
%     tem_logcoordinates=logm(orientation.Data(:,:,i)*R_Estimation(:,:,i)');
%     error_R(:,i)=skew2axis(tem_logcoordinates);
% end
% figure
% for i=1:3
%     subplot(3,1,i)
%     plot(2:N,error_R(i,2:N))
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
