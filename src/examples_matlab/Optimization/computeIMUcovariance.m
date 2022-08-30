function [Sigma_ik,Sigma_leg_ik] ...
    = computeIMUcovariance(i,j,w_meas,bar_bg,delta_bg,a_meas,bar_ba,delta_ba,Delta_t,sigma_g,simga_a,...
    sigma_legdrift,encoder_meas,contact_mat,encoder_dot_meas)
%computecovariance Compute covariance Sigma_ij for IMU and leg
% delta_bg ---- 


if j<i
    error("Wrong index for function computecovariance!")
elseif j==i
    Sigma_ik = zeros(9);
    warning("Index warning for function computecovariance!")
else
    Sigma_eta = blkdiag(sigma_g,simga_a);
    Sigma_ik = zeros(9);
    Sigma_leg_ik = zeros(3);
    Delta_R_ik = eye(3);
    A = eye(9);
    B = zeros(9,6);
    Mat_Cov_leg = zeros(3,9);
    for k = i:j-1

        % ceofficient
        A(1:3,1:3) = Exp((w_meas(:,k)-bar_bg-delta_bg)*Delta_t);
        A(4:6,1:3) = -Delta_R_ik * axis2skew(a_meas(:,k)-bar_ba-delta_ba)*Delta_t;
        A(7:9,1:3) = -0.5 * Delta_R_ik * axis2skew(a_meas(:,k)-bar_ba-delta_ba)*Delta_t*Delta_t;
        A(7:9,4:6) = Delta_t*eye(3);

        B(1:3,1:3) = J_R((w_meas(:,k)-bar_bg-delta_bg)*Delta_t) * Delta_t;
        B(4:6,4:6) = Delta_R_ik * Delta_t;
        B(7:9,4:6) = 0.5 * Delta_R_ik * Delta_t * Delta_t;
        
        % Leg part------------------------------
        alpha_l = floor(contact_mat(k,1));
        alpha_r = floor(contact_mat(k,2));
        if alpha_l ==0 && alpha_r ==0
            error("No contact!")
        end
        %l_ik
        tilde_v_k_left = -axis2skew(w_meas(:,k)-bar_bg-delta_bg)* p_VectorNav_to_LeftToeBottom(encoder_meas(k,:))...
            -J_VectorNav_to_LeftToeBottom(encoder_meas(k,:))*encoder_dot_meas(:,k);
        tilde_v_k_right = -axis2skew(w_meas(:,k)-bar_bg-delta_bg)*p_VectorNav_to_RightToeBottom(encoder_meas(k,:))...
                -J_VectorNav_to_RightToeBottom(encoder_meas(k,:))*encoder_dot_meas(:,k);
        tilde_v_k  = (alpha_l *tilde_v_k_left + alpha_r * tilde_v_k_right)/(alpha_l+alpha_r);
        f_alpha_left = J_VectorNav_to_LeftToeBottom(encoder_meas(k,:));
        f_alpha_right = J_VectorNav_to_RightToeBottom(encoder_meas(k,:));
        f_alpha_total_wedge =  (alpha_l*axis2skew(f_alpha_left)+alpha_r*axis2skew(f_alpha_right))/(alpha_l+alpha_r);

        Mat_Cov_leg(:,1:3) = Delta_R_ik * axis2skew(tilde_v_k) *Delta_t;
        Mat_Cov_leg(:,4:6) = Delta_R_ik * f_alpha_total_wedge *Delta_t;
        Mat_Cov_leg(:,4:6) = eye(3) * Delta_t;

        Delta_R_ik = Delta_R_ik * A(1:3,1:3);   % update Delta_R_ik

        % covariance
        Sigma_leg_ik = Sigma_leg_ik + Mat_Cov_leg * blkdiag(Sigma_ik(1:3,1:3),sigma_g,sigma_legdrift) * Mat_Cov_leg';
        Sigma_ik = A * Sigma_ik * A' + B * Sigma_eta * B';
    end
end




end