X_hat_Alternating = zeros(6,6,N);
bias_hat_Alternating = zeros(6,N);

cov_w_X = blkdiag(0.001*eye(3),0.001*eye(3),0.001*eye(3),0.001*eye(3));
cov_bg_all = blkdiag(1e-5*eye(3),1e-5*eye(3));
cov_v_Y = 0.01*eye(3);

M_b = [Delta_t*eye(6);zeros(6)];
F = eye(12);F(7:9,4:6) = Delta_t * eye(3);
H=[zeros(3,6),-eye(3),eye(3)];
d_const = [0 0 0 0 1 -1]';
D = [-eye(12),zeros(12);zeros(6,24);F,-eye(12);zeros(3,12),H];
%%
maxiter = 1;
for i =1:N
    if i ==1 
        X_k_check = X_true(:,:,1);
        bg_km1_hat = zeros(3,1);
        ba_km1_hat = zeros(3,1);
        Pk = 0.001*eye(18);
        D0 =[-eye(18);H,zeros(3,6)];
        W0 = blkdiag(0.001*eye(18),X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
        z = (X_k_check*Y_meas(:,i) - d_const);
        u = [zeros(18,1);z(1:3)];
        sol = (D0'/W0*D0)\D0'/W0*u;
        X_hat_Alternating(:,:,i) = expm(Rn2liealgebra(sol(1:12)))*X_k_check;
        Pk_inv = (D0'/W0*D0);
        Pk_KF = Pk - Pk*[H,zeros(3,6)]'/...
            ([H,zeros(3,6)]*Pk*[H,zeros(3,6)]'+X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)')...
            * [H,zeros(3,6)] *Pk;
    else
        %init
%         bg_km1 = bg_km1_hat;
%         ba_km1 = ba_km1_hat;
        buffer_last = zeros(24+12,1);
        delta_bg_km1 = zeros(3,1);
        delta_ba_km1 = zeros(3,1);
        for k = 1:maxiter
            % fix bias
            X_k_check = InEKFprediction(X_hat_Alternating(:,:,i-1),Delta_t,w_meas(:,i-1),a_meas(:,i-1),...
                bg_km1_hat + delta_bg_km1,ba_km1_hat+delta_ba_km1);
%             W = blkdiag(inv(Pk_inv),cov_w_X,X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
            W = blkdiag(Pk_KF,cov_w_X,X_k_check(1:3,1:3)*cov_v_Y*X_k_check(1:3,1:3)');
           
            z = (X_k_check*Y_meas(:,i) - d_const);
            u = [zeros(12,1);delta_bg_km1;delta_ba_km1;zeros(12,1);z(1:3)];
%             [~,~,~]=Schur_compute((D'/W*D),D'/W*u);
            sol_all = (D'/W*D)\D'/W*u;
            X_k_next = expm(Rn2liealgebra(sol_all(13:24)))*X_k_check;

            % fix states
            z_forbias = [sol_all(1:12);zeros(6,1);sol_all(13:24);zeros(6,1)];
%             W_forbias = blkdiag(inv(Pk_inv),cov_bg_all,cov_w_X);
            W_forbias = blkdiag(Pk_KF,cov_bg_all,cov_w_X);
            D_forbias = [zeros(12,12);-eye(6),zeros(6);Adjoint(X_k_next)*M_b,zeros(12,6);zeros(6,6),-eye(6)];
            sol_bias = (D_forbias'/W_forbias*D_forbias)\D_forbias'/W_forbias*z_forbias;
            delta_bg_km1 =  sol_bias(1:3);
            delta_ba_km1 =  sol_bias(4:6);

            if abs(sum(sol_bias(7:12)))>1e-9
                error("Here!")
            end

            if k>2 && norm([sol_all;sol_bias]-buffer_last)<1e-6
%                 disp("alternating converge at "+num2str(k)+"-th iteration")
                break;
            end
            if k == maxiter
%                 disp("maximum iteration!")
%                 pause
            end
            buffer_last(1:24) = sol_all;
            buffer_last(25:end) = sol_bias;
        end
        % update bias hat
        bg_km1_hat = bg_km1_hat + delta_bg_km1;
        ba_km1_hat = ba_km1_hat + delta_ba_km1;
        % update estimation 
        X_hat_Alternating(:,:,i) = X_k_next;
        bias_hat_Alternating(1:3,i) = bg_km1_hat + delta_bg_km1;
        bias_hat_Alternating(4:6,i) = ba_km1_hat+delta_ba_km1;
        %update covirance
        D_cov = [eye(12),zeros(12,6+12+6);
                 zeros(6,12),eye(6),zeros(6,12+6);
                 -F,Adjoint(Liegroup_inverse(X_k_next)) * M_b,eye(12),zeros(12,6);
                 zeros(6,12),-eye(6),zeros(6,12),eye(6);
                 zeros(3,12+6),-H,zeros(3,6)];
        W_cov = blkdiag(inv(Pk_inv),cov_w_X,cov_bg_all,X_k_next(1:3,1:3)*cov_v_Y*X_k_next(1:3,1:3)');
        P_all = D_cov' /W_cov *D_cov; 
        [~,Pk_inv,~]=Schur_compute(P_all,zeros(size(P_all,1),1));

        F_KF = eye(12+6);
        F_KF(7:9,4:6) = Delta_t * eye(3);
        F_KF(1:12,13:end) = Adjoint(Liegroup_inverse(X_k_next)) * M_b;
        Pk = F_KF*Pk_KF*F_KF'+ blkdiag(cov_w_X,cov_bg_all);
        Pk_KF = Pk - Pk*[H,zeros(3,6)]'/...
            ([H,zeros(3,6)]*Pk*[H,zeros(3,6)]'+X_k_next(1:3,1:3)*cov_v_Y*X_k_next(1:3,1:3)')...
            * [H,zeros(3,6)] *Pk;
%         if abs(1-norm(Pk_inv*Pk_KF))>1e-6
%             error("Pk error")
%         end
    end
end


% position
figure
for i=1:3
    subplot(3,1,i)
    temp_a = reshape(X_hat_Alternating(1:3,5,:),[3,N]);
    plot(1:N,temp_a(i,:));
    hold on
    temp_b = reshape(X_true(1:3,5,:),[3,N]);
    plot(1:N,temp_b(i,:));
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
sgtitle("Alternating Estimation of linear position")

disp("Alternating InEKF --- ATE of position: " +num2str(norm(temp_b-temp_a)))

% bias
figure
for i=1:2
    for j=1:3
        subplot(3,2,i+2*(j-1))     
        if i==1
            plot(bias_hat_Alternating(j,:))
            hold on
            plot(bg_true(j,:))
            hold off 
            legend('estimation','groudtruth')
            title('bg')
        else
            plot(bias_hat_Alternating(3+j,:))
            hold on
            plot(ba_true(j,:))
            hold off 
            legend('estimation','groudtruth')
            title('ba')
        end
    end
end
sgtitle("Alternating: estimation of Bias")


