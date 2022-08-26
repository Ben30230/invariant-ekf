function [Bar_Delta_R,Bar_Delta_v,Bar_delta_p,Bar_delta_l,par_R_par_bg,par_v_par_bg,par_v_par_ba,par_p_par_bg,par_p_par_ba,par_l_par_bg]...
    = DeltabarRvpij(i,j,w_meas,bar_bg,a_meas,bar_ba,Delta_t,encoder_meas,contact_mat,encoder_dot_meas)
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
    Bar_delta_l = zeros(3,1);

    % par_R_par_bg = zeros(3);
    par_v_par_bg = zeros(3);
    par_v_par_ba = zeros(3);
    par_p_par_bg = zeros(3);
    par_p_par_ba = zeros(3);
    par_l_par_bg = zeros(3);
    
    par_R_par_bg_ik = zeros(3);
    for k=i:j-1
        % Leg part------------------------------
        alpha_l = floor(contact_mat(k,1));
        alpha_r = floor(contact_mat(k,2));
        if alpha_l ==0 && alpha_r ==0
            error("No contact!")
        end
        %l_ik
        bar_v_k_left = -axis2skew(w_meas(:,k)-bar_bg)* p_VectorNav_to_LeftToeBottom(encoder_meas(k,:))...
            -J_VectorNav_to_LeftToeBottom(encoder_meas(k,:))*encoder_dot_meas(:,k);
        bar_v_k_right = -axis2skew(w_meas(:,k)-bar_bg)*p_VectorNav_to_RightToeBottom(encoder_meas(k,:))...
                -J_VectorNav_to_RightToeBottom(encoder_meas(k,:))*encoder_dot_meas(:,k);
        bar_v_k  = (alpha_l *bar_v_k_left + alpha_r * bar_v_k_right)/(alpha_l+alpha_r);
        Bar_delta_l = Bar_delta_l + Bar_Delta_R * bar_v_k * Delta_t;
        % par_l_par_bg
        f_alpha_left = J_VectorNav_to_LeftToeBottom(encoder_meas(k,:));
        f_alpha_right = J_VectorNav_to_RightToeBottom(encoder_meas(k,:));
        f_alpha_total_wedge =  (alpha_l*axis2skew(f_alpha_left)+alpha_r*axis2skew(f_alpha_right))/(alpha_l+alpha_r);
        par_l_par_bg = par_l_par_bg - Bar_Delta_R * (f_alpha_total_wedge+axis2skew(bar_v_k)*par_R_par_bg_ik) * Delta_t;
        %---------------------------------------------

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

