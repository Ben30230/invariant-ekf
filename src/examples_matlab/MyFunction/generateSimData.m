clc,clear
close all
load data\measurements\angular_velocity.mat
load data\measurements\linear_acceleration.mat
% generate simulation data

flag_bias = 1;

% input 
R0 = diag([1,-1,-1]);
v0 = zeros(3,1);
p0 = zeros(3,1);
d0 = zeros(3,1);

N=10000;
w_true = angular_velocity.Data(2:2:end,:)';
% a_true = [linear_acceleration.Data(200:396,:)',linear_acceleration.Data(396:2:end,:)'];
a_true = [linear_acceleration.Data(2:2:end,:)'];
%
cov_angular_vel = 0.00001*eye(3);
cov_acc = 0.00002*eye(3);
cov_bg = 1e-9*eye(3);
cov_ba = 1e-9*eye(3);
cov_yk = 0.000001*eye(3);

gravity = [0 0 -9.81]';
Delta_t = 0.001;

% output 
X_true = zeros(6,6,N+1);
w_meas = zeros(3,N);
a_meas = zeros(3,N);
bg_true = zeros(3,N+1);
ba_true = zeros(3,N+1);
Y_meas = zeros(6,N+1);

% Dynamcis 
X_true(:,:,1) = [R0,v0,p0,d0;zeros(3),eye(3)];
for i=1:N
    Gamma_i = [eye(3),Delta_t*gravity,0.5*gravity*Delta_t*Delta_t,zeros(3,1);zeros(3),eye(3)];
    Phi_i = X_true(:,:,i);
    Phi_i(1:3,5) = Phi_i(1:3,5)+Delta_t * Phi_i(1:3,4);
    Upsilon_i = [Exp(w_true(:,i)*Delta_t),a_true(:,i)*Delta_t,zeros(3,2);zeros(3),eye(3)];
    X_true(:,:,i+1) = Gamma_i*Phi_i*Upsilon_i;
    if flag_bias == 1
%         bg_true(:,i+1) = bg_true(:,i) + multirandn(zeros(3,1),cov_bg);
%         ba_true(:,i+1) = ba_true(:,i) + multirandn(zeros(3,1),cov_ba);
        bg_true(:,i+1) = 0.0005/(1+exp(-i*Delta_t));
        ba_true(:,i+1) = 0.0002/(1+exp(-i*Delta_t));
    end
end

%measurements
for i = 1:N+1
    if i < (N+1)
        w_meas(:,i) = w_true(:,i) + bg_true(:,i) + multirandn(zeros(3,1),cov_angular_vel);
        a_meas(:,i) = a_true(:,i) + ba_true(:,i) + multirandn(zeros(3,1),cov_acc);
    end
    Y_meas(:,i) = X_true(:,:,i)\[zeros(1,3) 0 1 -1]' + [multirandn(zeros(3,1),cov_yk);zeros(3,1)];
end

save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\w_true.mat',"w_true")
save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\a_true.mat',"a_true")
save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\X_true.mat',"X_true")
save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\bg_true.mat',"bg_true")
save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\ba_true.mat',"ba_true")

save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\w_meas.mat',"w_meas")
save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\a_meas.mat',"a_meas")
save('C:\Users\Ben\Desktop\Ph.D\Project\InEKF\invariant-ekf\src\examples_matlab\data\Y_meas.mat',"Y_meas")


%%
%plot 
% figure
% for i=1:3
%     subplot(3,1,i)
%     temp_a = reshape(X_true(1:3,4,:),[3,N+1]);
%     plot(temp_a(i,:));
%     legend('GroundTruth')
%     xlabel('time step')
%     switch i
%         case 1
%             ylabel('x(m)')
%         case 2
%             ylabel('y(m)')
%         case 3
%             ylabel('z(m)')
%     end
% end
% sgtitle("Estimation of linear position")

% plot(bg_true(1,:))



