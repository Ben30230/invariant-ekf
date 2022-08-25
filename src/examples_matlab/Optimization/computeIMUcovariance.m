function [Sigma_ik] = computeIMUcovariance(i,j,w_meas,bar_bg,delta_bg,a_meas,bar_ba,delta_ba,Delta_t,sigma_g,simga_a)
%computecovariance Compute covariance Sigma_ij for IMU
% delta_bg ---- 


if j<i
    error("Wrong index for function computecovariance!")
elseif j==i
    Sigma_ik = zeros(9);
    warning("Index warning for function computecovariance!")
else
    Sigma_eta = blkdiag(sigma_g,simga_a);
    Sigma_ik = zeros(9);
    Delta_R_ik = eye(3);
    A = eye(9);
    B = zeros(9,6);
    for k = i:j-1
        % ceofficient
        A(1:3,1:3) = Exp((w_meas(:,k)-bar_bg-delta_bg)*Delta_t);
        A(4:6,1:3) = -Delta_R_ik * axis2skew(a_meas(:,k)-bar_ba-delta_ba)*Delta_t;
        A(7:9,1:3) = -0.5 * Delta_R_ik * axis2skew(a_meas(:,k)-bar_ba-delta_ba)*Delta_t*Delta_t;
        A(7:9,4:6) = Delta_t*eye(3);

        B(1:3,1:3) = J_R((w_meas(:,k)-bar_bg-delta_bg)*Delta_t) * Delta_t;
        B(4:6,4:6) = Delta_R_ik * Delta_t;
        B(7:9,4:6) = 0.5 * Delta_R_ik * Delta_t * Delta_t;

        Delta_R_ik = Delta_R_ik * A(1:3,1:3);   % update Delta_R_ik

        % covariance
        Sigma_ik = A * Sigma_ik * A' + B * Sigma_eta * B';
    end
end




end