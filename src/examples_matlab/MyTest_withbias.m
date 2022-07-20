clc,clear
close all
addpath('MyFunction\')
addpath(genpath('forward_kinematics'))

load data\measurements\angular_velocity.mat
load data\measurements\linear_acceleration.mat
load data\measurements\contact.mat
load data\measurements\encoders.mat
load data\measurements\X_init.mat

t = angular_velocity.Time;
angular_mat = angular_velocity.Data;
acc_mat = linear_acceleration.Data;
joint_meas=encoders.Data;
contact_mat = contact.Data;

% ture state
load data\ground_truth\orientation.mat
load data\ground_truth\position.mat
load data\ground_truth\velocity.mat

%%
clc,clear
g=[0 0 -9.81];
v=rand(3,1);
p=rand(3,1);
dl=rand(3,1);
dr=rand(3,1);
R = expm(axis2skew(rand(3,1)));

C=rand(21);
P=C'*C;

At1=[zeros(3,15),-R,zeros(3);
                        axis2skew(g),zeros(3,12),-axis2skew(v)*R,-R;
                        zeros(3),eye(3),zeros(3,9),-axis2skew(p)*R,zeros(3);
                        zeros(3,15),-axis2skew(dl)*R,zeros(3);
                        zeros(3,15),-axis2skew(dr)*R,zeros(3);
                        zeros(6,21)];            %21*21


At1*P*At1'

At2 = [zeros(3,12),-R,zeros(3);
                        axis2skew(g),zeros(3,9),-axis2skew(v)*R,-R;
                        zeros(3),eye(3),zeros(3,6),-axis2skew(p)*R,zeros(3);
                        zeros(3,12),-axis2skew(dl)*R,zeros(3);
                        zeros(3,12),-axis2skew(dr)*R,zeros(3);
                        zeros(6,18)];            %18






