function [Delta_l_ij] = Deltatilde_l_ij(i,j,w_meas,bar_bg,delta_bg,Delta_t,tilde_v_k)
%Deltatilde_l_ij For leg residual function
%   为了简单，这里先分开写，其实可以和 IMU 部分整合到一起

Delta_l_ij = zeros(3,1);
Delta_R_ik = eye(3);
for k = i:j-1
%      = -axis2skew(w_meas(:,k)-bar_bg-delta_bg)*
    Delta_l_ij = Delta_l_ij + Delta_R_ik * tilde_v_k *Delta_t;
    Delta_R_ik = Delta_R_ik * Exp((w_meas(:,k)-bar_bg-delta_bg)*Delta_t);
end



end