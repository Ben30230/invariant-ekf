function [JacobianR,Jacobianv,Jacobianp] = JacobianofIMUresidualfunc(i,j,Ri,Rj,vi,vj,pi,pj,residual_R,delta_bg,Delta_t,...
    par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba,g)
%Jacobianofresidualfunc Compute jacobians of three residual functions
% input:
%       i ---- int
%       j ---- int
%       Ri ---- SO(3) current estimated Ri
%       Rj ---- SO(3) (current estimated)
%       vi ---- R^3 (current estimated)
%       vj ---- R^3 (current estimated)
%       pi ---- R^3 (current estimated)
%       pj ---- R^3 (current estimated)
%       residual_R ---- current residual function value of Delta orientation Rij
%       delta_bg ---- current value of bg
%       Delta_t ---- cycle time
%       par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba
%       ---- from function called "DeltabarRvpij"
%       g ---- gravity
% outout:
%       JacobianR: 3*24 mat order: phii vi pi bgi bai phij vj pj
%       Jacobianv: 3*24 mat order: phii vi pi bgi bai phij vj pj
%       Jacobianp: 3*24 mat order: phii vi pi bgi bai phij vj pj



if j<i
    error("Wrong index for function computecovariance!")
elseif j==i
    JacobianR = zeros(3,24);
    Jacobianv = zeros(3,24);
    Jacobianp = zeros(3,24);
    warning("Index warning for function computecovariance!")
else
    Deltatij = (j-i) * Delta_t;
    Deltatij_square = (j-i) * (j-i) * Delta_t * Delta_t;

    % position residual function jacobian
    par_rp_par_phii = axis2skew(Ri'*(pj-pi-vi*Deltatij-0.5*9.81*Deltatij_square));
    par_rp_par_vi = -Ri'*Deltatij;
    par_rp_par_pi = -eye(3);
    par_rp_par_bg = -par_p_par_bg;
    par_rp_par_ba = -par_p_par_ba;
    par_rp_par_pj = Ri'*Rj;

    Jacobianp = [par_rp_par_phii,par_rp_par_vi,par_rp_par_pi,par_rp_par_bg,...
        par_rp_par_ba,zeros(3,6),par_rp_par_pj];

    % velocity residual function jacobian
    par_rv_par_phii = axis2skew(Ri'*(vj-vi-g*Deltatij));
    par_rv_par_vi = -Ri';
    par_rv_par_bg = -par_v_par_bg;
    par_rv_par_ba = -par_v_par_ba;
    par_rv_par_vj = Ri';
    
    Jacobianv = [par_rv_par_phii,par_rv_par_vi,zeros(3),par_rv_par_bg,...
        par_rv_par_ba,zeros(3),par_rv_par_vj,zeros(3)];

    % orientation residual function jacobian
    par_rR_par_phii = -J_R_inv(residual_R)*Rj'*Ri;
    par_rR_par_bg = -J_R_inv(residual_R) * Exp(residual_R)'*J_R(par_R_par_bg*delta_bg)*par_R_par_bg;
    par_rR_par_phij = J_R_inv(residual_R);
    
    JacobianR = [par_rR_par_phii,zeros(3,6),par_rR_par_bg,...
        zeros(3),par_rR_par_phij,zeros(3,6)];
end



end