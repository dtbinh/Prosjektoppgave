
N_measurment                = 8;

% Uncertainties               
if cfg.M_Uncertainty == 1
    r_PR            = 0.01;
else
    r_PR            = 0;
end
R_PR                = eye(N_measurment) * r_PR;

q_PR                = 1;
Q_PR                = 1*q_PR*diag([0,0,0,10,10,2,0.01,0.1,0.1]);

%System

T                   = 100000;                %Forgetting factor
% cfg.A_EKF           = [zeros(3,3),      eye(3) ,  zeros(3,3);
%                        zeros(3,3), eye(3)*-d_t ,     -eye(3);
%                        zeros(3,3),  zeros(3,3) , -1/T*eye(3)];
                   
cfg.A_EKF           = [zeros(3,3),      eye(3) ,  zeros(3,3);
                       zeros(3,3), eye(3)*-d_t ,     -eye(3);
                       zeros(3,3),  zeros(3,3) ,  zeros(3,3)];
                   
cfg.B_EKF           = [zeros(3,3);
                           eye(3);
                       zeros(3,3)];

%Bias
cfg.bias_am         = 1 * rand(3,1);
log.add('bias_am',3);
log.store('bias_am', bias_m, 1);

[A_EKF, B_EKF]      = Convert_cont2disc(cfg.A_EKF, cfg.B_EKF, h);

%Radio Beacons
sig_pos_1                   = [ 1;  1;  1];
sig_pos_2                   = [ 1;  1; -1];
sig_pos_3                   = [ 1; -1;  1];
sig_pos_4                   = [-1;  1;  1];
sig_pos_5                   = [ 1; -1; -1];
sig_pos_6                   = [-1; -1;  1];
sig_pos_7                   = [-1;  1; -1];
sig_pos_8                   = [-1; -1; -1];
sig_pos                     = 30*[sig_pos_1, sig_pos_2, sig_pos_3, sig_pos_4, sig_pos_5, sig_pos_6, sig_pos_7, sig_pos_8];
clear sig_pos_1 sig_pos_2 sig_pos_3 sig_pos_4 sig_pos_5 sig_pos_6 sig_pos_7 sig_pos_8

%H                           = diag([1, 1, 1, 0, 0, 0, 0, 0, 0]);

H_EK                        = @(x, sig_pos, length) [(x-sig_pos(:,1))'/length(1), zeros(1,3), zeros(1,3);
                                                     (x-sig_pos(:,2))'/length(2), zeros(1,3), zeros(1,3);
                                                     (x-sig_pos(:,3))'/length(3), zeros(1,3), zeros(1,3);
                                                     (x-sig_pos(:,4))'/length(4), zeros(1,3), zeros(1,3);
                                                     (x-sig_pos(:,5))'/length(5), zeros(1,3), zeros(1,3);
                                                     (x-sig_pos(:,6))'/length(6), zeros(1,3), zeros(1,3);
                                                     (x-sig_pos(:,7))'/length(7), zeros(1,3), zeros(1,3);
                                                     (x-sig_pos(:,8))'/length(8), zeros(1,3), zeros(1,3)];
