clc
clear all
close all

%%

p_e_eb1 = [500000; 500000; -500000];
p_e_eb2 = [500000; -500000; -500000];
p_e_eb3 = [-500000; 500000; -500000];
p_e_eb4 = [-500000; -500000; -500000];
p_e_eb5 = [0; -1000000; -500000];
p_e_eb6 = [0; 1000000; -500000];
p_e_eb7 = [-1000000; 0; -500000];
p_e_eb8 = [1000000; 0; -500000];

p_e_eb11 = [100; 100; -10];
p_e_eb12 = [-100; 100; -3];
p_e_eb13 = [-100; 100; -10];
p_e_eb14 = [-100; -100; -7];

% true clock error
beta_gnss = 4;
beta_uwb = 7;

% true pos
p_e_eb = [20; -30; -50];

% intit estimate
p_e_eb0 = [0; 0; 0];
s_states.pos = p_e_eb0; 
s_states.velo = [0;0;0];
s_states.beta_gnss = 0;
s_states.beta_uwb = 0;
s_states.R_en = eye(3);
% P_ = blkdiag( 10*eye(3), 0.1*eye(3), 10, 10);
P_ = blkdiag( 10*eye(3), 0.1*eye(3), 10);

% meas noise 
std_ = 1;
std_scale = 1;

% process noise
g = 9.81;
std_acc = 0.1*g;
std_beta_gnss = 0.01;
std_beta_uwb = 0.05;
% Q = blkdiag( std_acc^2*eye(3), std_beta_gnss^2, std_beta_uwb^2);
Q = blkdiag( std_acc^2*eye(3), std_beta_gnss^2);

c_gnss_meas = cell(0);
c_uwb_meas = cell(0);
%% 1
s_gnss.p_eb_e_k = p_e_eb1;
rho_gnss(1) = norm(p_e_eb-p_e_eb1);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%% 2
s_gnss.p_eb_e_k = p_e_eb2;
rho_gnss(2) = norm(p_e_eb-p_e_eb2);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%% 3
s_gnss.p_eb_e_k = p_e_eb3;
rho_gnss(3) = norm(p_e_eb-p_e_eb3);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%% 4
s_gnss.p_eb_e_k = p_e_eb4;
rho_gnss(4) = norm(p_e_eb-p_e_eb4);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%% 5
s_gnss.p_eb_e_k = p_e_eb5;
rho_gnss(5) = norm(p_e_eb-p_e_eb5);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%% 6
s_gnss.p_eb_e_k = p_e_eb6;
rho_gnss(6) = norm(p_e_eb-p_e_eb6);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%% 7
s_gnss.p_eb_e_k = p_e_eb7;
rho_gnss(7) = norm(p_e_eb-p_e_eb7);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%% 8
s_gnss.p_eb_e_k = p_e_eb8;
rho_gnss(8) = norm(p_e_eb-p_e_eb8);
s_gnss.std = std_*std_scale;
c_gnss_meas{length(c_gnss_meas)+1} = s_gnss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% 1
% s_uwb.p_eb_e_k = p_e_eb11;
% rho_uwb(1) = norm(p_e_eb-p_e_eb11);
% s_uwb.std = std_*std_scale;
% c_uwb_meas{length(c_uwb_meas)+1} = s_uwb;
% 
% %% 2
% s_uwb.p_eb_e_k = p_e_eb12;
% rho_uwb(2) = norm(p_e_eb-p_e_eb12);
% s_uwb.std = std_*std_scale;
% c_uwb_meas{length(c_uwb_meas)+1} = s_uwb;
% 
% %% 3
% s_uwb.p_eb_e_k = p_e_eb13;
% rho_uwb(3) = norm(p_e_eb-p_e_eb13);
% s_uwb.std = std_*std_scale;
% c_uwb_meas{length(c_uwb_meas)+1} = s_uwb;
% 
% %% 4
% s_uwb.p_eb_e_k = p_e_eb14;
% rho_uwb(4) = norm(p_e_eb-p_e_eb14);
% s_uwb.std = std_*std_scale;
% c_uwb_meas{length(c_uwb_meas)+1} = s_uwb;

N_meas_gnss = length(c_gnss_meas);
N_meas_uwb = length(c_uwb_meas);
C0_bench = [];
for k3 = 1:N_meas_gnss
    C0_bench = [C0_bench; (p_e_eb - c_gnss_meas{k3}.p_eb_e_k)'/(rho_gnss(k3)) 1 0];
end
for k4 = 1:N_meas_uwb
    C0_bench = [C0_bench; (p_e_eb - c_uwb_meas{k4}.p_eb_e_k)'/(rho_uwb(k4)) 0 1]; 
end

if N_meas_uwb == 0
    C0_bench(:,end) = [];
end
if N_meas_gnss == 0
    C0_bench(:,4) = [];
end
G_bench = inv(C0_bench'*C0_bench);
HDOP_bench = sqrt( G_bench(1,1) + G_bench(2,2) )
VDOP_bench = sqrt( G_bench(3,3) )
TDOP_GNSS_bench = sqrt( G_bench(4,4) )
if N_meas_uwb > 0
    TDOP_UWB_bench = sqrt( G_bench(5,5) )
end


%%
nav_frame = 'ned';
Ts = 0.2;
num_of_sec_of_sim = 1000;
N = num_of_sec_of_sim/Ts;

pos_hat_data = zeros(3,N);
beta_gnss_hat_data = zeros(1,N);
beta_uwb_hat_data = zeros(1,N);
pos_data = p_e_eb.*ones(3,N);

hdop = HDOP_bench*ones(1,N);
hdop_hat = zeros(1,N);
vdop = VDOP_bench*ones(1,N);
vdop_hat = zeros(1,N);
tdop_gnss = TDOP_GNSS_bench*ones(1,N);
tdop_gnss_hat = zeros(1,N);

if N_meas_uwb > 0
    tdop_uwb = TDOP_UWB_bench*ones(1,N);
else
    tdop_uwb = zeros(1,N);
end
tdop_uwb_hat = zeros(1,N);

beta_gnss_data = beta_gnss*ones(1,N);
beta_uwb_data = beta_uwb*ones(1,N);

% noise gen
noise_gnss = normrnd( 0, std_, N_meas_gnss, N);
noise_uwb = normrnd( 0, std_, N_meas_uwb, N);
y_tilde_data =  zeros(N_meas_gnss, N);
% y_tilde_data =  zeros(N_meas_gnss + N_meas_uwb, N);
for k1 = 1:N
    for k2 = 1:N_meas_gnss
        c_gnss_meas{k2}.y = rho_gnss(k2) + noise_gnss( k2, k1) + beta_gnss;
    end
    
    for k3 = 1:N_meas_uwb
        c_uwb_meas{k3}.y = rho_uwb(k3) + noise_uwb( k3, k1) + beta_uwb; 
    end
    [ s_states, P_, y_tilde ] = step_gnss( Ts, s_states, P_, c_gnss_meas, Q, nav_frame );
%     [ s_states, P_, y_tilde ] = step_gnss_uwb( Ts, s_states, P_, c_gnss_meas, c_uwb_meas, Q, nav_frame );
    pos_hat_data(:,k1) = s_states.pos;
    beta_gnss_hat_data(k1) = s_states.beta_gnss;
%     beta_uwb_hat_data(k1) = s_states.beta_uwb;
    hdop_hat(k1) = s_states.HDOP;
    vdop_hat(k1) = s_states.VDOP;
    tdop_gnss_hat(k1) = s_states.TDOP_GNSS;
%     tdop_uwb_hat(k1) = s_states.TDOP_UWB;
    
    % if size of y_tilde changes duo to satelittes with to low elevation a more elegant implemetation is needed
    % this is not an issue in the NED simulation
    y_tilde_data(:,k1) = y_tilde; 
end

c_leg = {'True', 'Est'};
figure
subplot(3,1,1)
    plot(pos_data(1,:))
    hold on;
    plot(pos_hat_data(1,:),'--')
    hold off;
    grid on;
    ylabel('x [m]')
    title('Position and position estimate')
    legend(c_leg,'location','best')
subplot(3,1,2)
    plot(pos_data(2,:))
    hold on;
    plot(pos_hat_data(2,:),'--')
    hold off;
    grid on;
    ylabel('y [m]')
subplot(3,1,3)
    plot(pos_data(3,:))
    hold on;
    plot(pos_hat_data(3,:),'--')
    hold off;
    grid on;
    ylabel('z [m]')
    xlabel('Time [sec]')

c_leg = {'beta gnss', 'est beta gnss'};
figure
plot(beta_gnss_data)
hold on;
plot(beta_gnss_hat_data)
title('clock error estimate')
grid on;
legend(c_leg, 'location','best')

c_leg = {'beta uwb', 'est beta uwb'};
figure
plot(beta_uwb_data)
hold on;
plot(beta_uwb_hat_data)
title('clock error estimate')
grid on;
legend(c_leg, 'location','best')
 
c_leg = {'True DOP', 'Est DOP'};
figure
subplot(4,1,1)
    plot(hdop)
    hold on;
    plot(hdop_hat, '--')
    hold off;
    grid on;
    title('HDOP')
    legend(c_leg,'location','best')
subplot(4,1,2)
    plot(vdop)
    hold on;
    plot(vdop_hat, '--')
    hold off;
    grid on;
    title('VDOP')
subplot(4,1,3)
    plot(tdop_gnss)
    hold on;
    plot(tdop_gnss_hat, '--')
    hold off;
    grid on;
    title('TDOP GNSS')
subplot(4,1,4)
    plot(tdop_uwb)
    hold on;
    plot(tdop_uwb_hat, '--')
    hold off;
    grid on;
    title('TDOP UWB')
    
if (N_meas_gnss > 0) && (N_meas_uwb > 0)
    figure
    subplot(2,1,1)
    hold on;
    for k = 1:N_meas_gnss
        histogram( y_tilde_data(k,:) )
    end
    hold off;
    grid on;
    title('y-tilde histogram gnss')

    subplot(2,1,2)
    hold on;
    for k = N_meas_gnss+1:N_meas_gnss+N_meas_uwb
        histogram( y_tilde_data(k,:) )
    end
    hold off;
    grid on;
    title('y-tilde histogram uwb')
elseif N_meas_gnss > 0
    figure
    hold on;
    for k = 1:N_meas_gnss
        histogram( y_tilde_data(k,:) )
    end
    hold off;
    grid on;
    title('y-tilde histogram gnss') 
elseif N_meas_uwb > 0
    figure
    hold on;
    for k = 1:N_meas_gnss
        histogram( y_tilde_data(k,:) )
    end
    hold off;
    grid on;
    title('y-tilde histogram uwb') 
end
    
rms_x = rms(pos_data(1, ceil(end/2):end )-pos_hat_data(1, ceil(end/2):end ) )
rms_y = rms(pos_data(2, ceil(end/2):end )-pos_hat_data(2, ceil(end/2):end ) )
rms_z = rms(pos_data(3, ceil(end/2):end )-pos_hat_data(3, ceil(end/2):end ) )