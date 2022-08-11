clc,clear  
g = [0 0 -9.81]';
At = [zeros(3,15)
    axis2skew(g),zeros(3,12)
    zeros(3),eye(3),zeros(3,9)
    zeros(3,12),eye(3)
    zeros(3,15)];
At*At
