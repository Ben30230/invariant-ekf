function [output] = calculatepar_R_par_bgij(i,j,w_meas,bar_bg,Delta_t)
%calculatepar_R_par_bgij

par_R_par_bg = zeros(3);

DealtaR_kplus1j = eye(3);
if (j-1)<i

else
    for k = (j-1):-1:i
        par_R_par_bg = par_R_par_bg - DealtaR_kplus1j' * J_R((w_meas(:,k)-bar_bg)*Delta_t)*Delta_t;
        DealtaR_kplus1j = Exp((w_meas(:,k)-bar_bg)*Delta_t) * DealtaR_kplus1j;
    end
end

output = par_R_par_bg;




end

