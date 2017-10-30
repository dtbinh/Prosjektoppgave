clear;
clf;
%% Switches

cfg.Gravity         = 1;
cfg.P_Uncertainty   = 1;
cfg.M_Uncertainty   = 1;
cfg.Integral        = 0;
cfg.SaveLog         = 0;
cfg.Model_3D        = 0;

%% Parameters

% Simulation parameters
h                   = 0.02;
tend                = 200;
time                = 0:h:tend-h;

% Observer

% Physical
if cfg.Gravity == 1
    g_z             = 9.81;
else
    g_z             = 0;
end
cfg.g = [0;0;g_z];
cfg.m_c             = 1;
m_c                 = cfg.m_c;

d_t                 = 0.1;
d_r                 = 0.1;

% System

cfg.A_Ctrl          = [zeros(3,3),      eye(3);
                       zeros(3,3), eye(3)*-d_t];

cfg.B_Ctrl          = [0, 0, 0;
                       0, 0, 0;
                       0, 0, 0;
                       1, 0, 0;
                       0, 1, 0;
                       0, 0, 1];

cfg.C               = [1, 0, 0, 0, 0, 0;
                       0, 1, 0, 0, 0, 0;
                       0, 0, 1, 0, 0, 0];
              
bias_m              = randn(3,1);
a                   = zeros(3,1);
if cfg.Gravity == 1
    cfg.k_f         = 4 * cfg.m_c * g_z;
else
    cfg.k_f         = 4 * cfg.m_c * 9.81;
end

cfg.max_ang         = pi/3;

%References

p_d             = [cos(time/12)-1; sin(time/12); cos(time/9)+time/24-1];
%p_d             = zeros(3,size(time,2));
%p_d             = [ones(size(time)); exp(-time/5); zeros(size(time))];
%p_d             = [0*ones(3,size(time,2)/3), 1*ones(3,size(time,2)/3), 5*ones(3,size(time,2)/3)];

T_p_d           = 0.05;
calc_diff       = @ (p_d_post, p_d_pre, h) (p_d_post - p_d_pre)/h; 

%% Variables
% References

log = Logger();
log.init(time,h);

%DOF variables and their derivates
log.add('p',3);
log.add('p_m',3);
log.add('v',3);
log.add('theta',3);

%Regulator Variables
log.add('tau',3);
log.add('M',3);

%References
log.add('p_d',3);
log.add('v_d',3);

log.add('eta_d',3);

%Measurements
log.add('x_hat', 3);
log.add('x_hat_filtered', 3);
log.add('error',8);
