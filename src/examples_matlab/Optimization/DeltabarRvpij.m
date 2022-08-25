function [Bar_Delta_R,Bar_Delta_v,Bar_delta_p,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba]...
    = DeltabarRvpij(i,j,w_meas,bar_bg,a_meas,bar_ba,Delta_t)
%DeltatildeRvpij  \prod_{k=i}^{j-1}\Exp\left( (\tilde{w}_k-b_i^g)\Delta t \right) if  j>i\\
% input ---- i int
%       ---- j int
%       ---- w_meas 3*(j-i)-dim mat
%       ---- bg 3-dim vec
%       ---- a_meas 3*(j-1)-dim mat
%       ---- ba 3-dim vec
%       ---- Delta_t scaler double
%
% output ---- SO(3)

if j < i
    error("Wrong Index of Delta Rvp!")
elseif j==i
    Bar_Delta_R = eye(3);
    Bar_Delta_v = zeros(3,1);
    Bar_delta_p = zeros(3,1);
    warning("i==j!")
else
    Bar_Delta_R = eye(3);
    Bar_Delta_v = zeros(3,1);
    Bar_delta_p = zeros(3,1);
%     par_R_par_bg = zeros(3);
    par_v_par_bg = zeros(3);
    par_v_par_ba = zeros(3);
    par_p_par_bg = zeros(3);
    par_p_par_ba = zeros(3);

    
%     par_R_par_bg=calculatepar_R_par_bgij(i,j,w_meas,bar_bg,Delta_t);
    par_R_par_bg_ik = zeros(3);
    for k=i:j-1
        %par_p_par_bg_ik
        par_p_par_bg = par_p_par_bg +par_v_par_bg*Delta_t-0.5*Bar_Delta_R*axis2skew(a_meas(:,k)-bar_ba)...
            *par_R_par_bg_ik*Delta_t^2;  
        %par_p_par_ba_ik
        par_p_par_ba = par_p_par_ba + par_v_par_ba*Delta_t - 0.5*Bar_Delta_R*Delta_t^2;
        %par_v_par_bg_ik
        par_v_par_bg = par_v_par_bg - Bar_Delta_R * axis2skew(a_meas(:,k)-bar_ba)...
            *par_R_par_bg_ik*Delta_t;
        %par_v_par_ba_ik
        par_v_par_ba = par_v_par_ba - Bar_Delta_R *Delta_t;

%         if norm(calculatepar_R_par_bgij(i,k,w_meas,bar_bg,Delta_t)-par_R_par_bg_ik,"fro")>1e-6
%             error("!")
%         end

        % update par_R_par_bg_ik
        par_R_par_bg_ik = Exp((w_meas(:,k)-bar_bg)*Delta_t)'*par_R_par_bg_ik - J_R((w_meas(:,k)-bar_bg)*Delta_t)*Delta_t;

        %p_ik
        Bar_delta_p = Bar_delta_p + Bar_Delta_v*Delta_t+0.5*Bar_Delta_R*(a_meas(:,k)-bar_ba)*Delta_t^2;
        %v_ik
        Bar_Delta_v = Bar_Delta_v + Bar_Delta_R*(a_meas(:,k)-bar_ba)*Delta_t;
        %R_ik
        Bar_Delta_R = Bar_Delta_R*Exp((w_meas(:,k)-bar_bg)*Delta_t);
    end
    % compute par_R_par_bg ij
    par_R_par_bg = par_R_par_bg_ik;
end

end

