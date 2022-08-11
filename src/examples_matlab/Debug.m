clc,clear
addpath('MyFunction\')
addpath(genpath('forward_kinematics'))
load data\measurements\angular_velocity.mat
load data\measurements\linear_acceleration.mat
load data\measurements\encoders.mat

load data\ground_truth\orientation.mat
load data\ground_truth\position.mat
load data\ground_truth\velocity.mat

%% Ture states 
angular_mat = angular_velocity.Data;
acc_mat = linear_acceleration.Data;
joint_meas=encoders.Data;

w_meas = angular_mat(9202:9526,:)';
acc_meas = acc_mat(9202:9526,:)';

bg_true = [0.03,0.03,0.03]';
ba_true = [0.02,0.02,0.02]';

w_true = w_meas-bg_true-multirandn(zeros(3,1),0.1*eye(3));
acc_true = acc_meas - ba_true - multirandn(zeros(3,1),0.1*eye(3));

N=length(w_true);

R_true = zeros(3,3,N);
v_true = zeros(3,N);
p_true = zeros(3,N);
d_true = zeros(3,N);
y_true = zeros(3,N);

alpha = joint_meas(9202:9526,:);
R_true(:,:,1) = orientation.Data(:,:,9202);
v_true(:,1) = velocity.Data(9202,:)';
p_true(:,1) = position.Data(9202,:);
d_true(:,1) = p_true(:,1) + R_true(:,:,1) * p_VectorNav_to_LeftToeBottom (alpha(1,:));


Delta_T = 0.0005;
g=[0,0,-9.81]';
for i = 2:N
    R_true(:,:,i) = R_true(:,:,i-1) * expm(axis2skew(w_true(:,i))*Delta_T);
    v_true(:,i) = v_true(:,i-1) + (R_true(:,:,i) * acc_true(:,i) + g )*Delta_T;
    p_true(:,i) = p_true(:,i-1) + v_true(:,i) * Delta_T + 0.5 * (R_true(:,:,i) * acc_true(:,i)+g)*Delta_T^2;
    d_true(:,i) = d_true(:,i-1) +multirandn(zeros(3,1),0.1*eye(3))*Delta_T ;

    y_true(:,i) = R_true(:,:,i)' * (d_true(:,i)-p_true(:,i)) + multirandn(zeros(3,1),0.01*eye(3));
end

%% InEKF
R_hat = zeros(3,3,N);
v_hat = zeros(3,N);
p_hat = zeros(3,N);
d_hat = zeros(3,N);
bg_hat = zeros(3,N);
ba_hat = zeros(3,N);

R_hat(:,:,1) = R_true(:,:,1);
v_hat(:,1) = v_true(:,1);
p_hat(:,1) = p_true(:,1);
d_hat(:,1) = d_true(:,1);

P = 0.001*eye(18);

for i = 2:N
    % prediction
    R_hat(:,:,i) = R_hat(:,:,i-1) * expm(axis2skew(w_meas(:,i)-bg_hat(:,i-1))*Delta_T);
    v_hat(:,i) = v_hat(:,i-1) + (R_hat(:,:,i) * (acc_meas(:,i)-ba_hat(:,i-1)) + g )*Delta_T;
    p_hat(:,i) = p_hat(:,i-1) + v_hat(:,i) * Delta_T + 0.5 * (R_hat(:,:,i) * (acc_meas(:,i)-ba_hat(:,i-1)) + g )*Delta_T^2;
    d_hat(:,i) = d_hat(:,i-1);

    X = [R_hat(:,:,i),v_hat(:,i),p_hat(:,i),d_hat(:,i);zeros(3),eye(3)];

    At=[zeros(3,12), -R_hat(:,:,i),zeros(3);
        axis2skew(g),zeros(3,9),-axis2skew(v_hat(:,i))*R_hat(:,:,i), - R_hat(:,:,i);
        zeros(3),eye(3),zeros(3,6),-axis2skew(p_hat(:,i))*R_hat(:,:,i),zeros(3);
        zeros(3),eye(3),zeros(3,6),-axis2skew(d_hat(:,i))*R_hat(:,:,i),zeros(3);
        zeros(6,18)];
    F=eye(size(At))+At*Delta_T;

    FQ = blkdiag(Adjoint(X),eye(6));
    h_R = R_VectorNav_to_LeftToeBottom(alpha(i,:));
    Q = blkdiag(0.01*eye(3),0.01*eye(3),zeros(3),h_R*0.01*eye(3)*h_R',0.01*eye(3),0.01*eye(3));
    P = F *P *F' + FQ * Q *FQ';

    % correction
    H = [zeros(3,6),-eye(3),eye(3),zeros(3,6)];
    R_prime = R_hat(:,:,i) * J_VectorNav_to_LeftToeBottom(alpha(i,:)) * 0.01*eye(14) ...
              * J_VectorNav_to_LeftToeBottom(alpha(i,:))' *R_hat(:,:,i)';
    K = P * H' /(H*P*H'+R_prime);
    delta = K * [eye(3),zeros(3)] * X * [y_true(:,i);0;1;-1];

    X=expm(Rn2liealgebra(delta(1:end-6)))*X;
    R_hat(:,:,i) = X(1:3,1:3);
    v_hat(:,i) = X(1:3,4);
    p_hat(:,i) = X(1:3,5);
    d_hat(:,i) = X(1:3,6);

    bg_hat(:,i) = bg_hat(:,i) + delta(end-5:end-3);
    ba_hat(:,i) = ba_hat(:,i) + delta(end-2:end);

    P = (eye(size(P))-K*H)*P;
    
end


%%
% position
figure
for i=1:3
    subplot(3,1,i)
    plot(1:N,p_hat(i,:))
    hold on
    plot(1:N,p_true(i,:))
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

% bias
figure
for i=1:2
    for j=1:3
        subplot(3,2,i+2*(j-1))     
        if i==1
            plot(1:N,bg_hat(j,:))
            title('bg')
        else
            plot(1:N,ba_hat(j,:))
            title('ba')
        end
    end
end
sgtitle("Estimation of Bias")







function AdjX = Adjoint(X)
            % Adjoint of SE_N(3)         
            N = size(X,2)-3;
            R = X(1:3,1:3);
            R_cell = repmat({R}, 1, N+1); 
            AdjX = blkdiag(R_cell{:});
            for i = 1:N
                AdjX(3*i+1:3*i+3,1:3) = axis2skew(X(1:3,i+3))*R;
            end
        end

function output = Rn2liealgebra(xi)
            xi=xi(:);
            N=length(xi)/3+2;
            output=zeros(N);
            for i=3:N
                if i==3
                    output(1:3,1:3)=axis2skew(xi(1:3));
                else
                    output(1:3,i)=xi(((i-3)*3+1):((i-3)*3+3));
                end
            end
        end



