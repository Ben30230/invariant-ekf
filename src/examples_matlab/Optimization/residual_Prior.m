function [r_Ri,r_vi,r_pi,r_bgi,r_bai] = residual_Prior(Ri,vi,pi,delta_bg,delta_ba,check_Ri,check_vi,check_pi,check_bgi,check_bai)
%residual_Prior
%
r_Ri = Log(check_Ri'*Ri);
r_vi = vi - check_vi;
r_pi = pi - check_pi;
r_bgi = delta_bg;
r_bai = delta_ba;

end