function [x, eta, eta_d, u, M, a] = Helicopter_3DOF(cfg, x, x_hat, eta, eta_d, p_d, v_d, p_integral, h)
%% Paramaeters
    
%Physical
m_c             = cfg.m_c;

d_t             = 0.1;
d_r             = 0.1;
g               = cfg.g;
i               = m_c;
I               = i * eye(3);

% System

A               = cfg.A_Ctrl;
               
B               = cfg.B_Ctrl;
               
[~,J1,J2]       = eulerang(eta(1), eta(2), eta(3));
               
%Regulators

%Transitional Regulator
omega_n_t       = 1.5;
zeta_t          = 1.2;

k_p_t           = m_c * omega_n_t^2;
k_d_t           = m_c * omega_n_t * 2 * zeta_t;
k_i_t           = m_c * omega_n_t * 0.1;

K_p_t           = k_p_t * eye(3);
K_d_t           = k_d_t * eye(3);
K_i_t           = k_i_t * eye(3);
K_t             = 1;

if cfg.Integral == 1
    tau             = @(x,p_integral) d_t * x(4:6) -K_t * (K_d_t * ( x(4:6) - v_d ) + K_p_t * ( x_hat(1:3) - p_d)) + K_i_t * (p_integral);
else
    tau             = @(x,p_integral) d_t * x(4:6) -K_t * (K_d_t * ( x(4:6) - v_d ) + K_p_t * ( x_hat(1:3) - p_d));
end

%Angular Regulator
omega_n_r       = 20;
zeta_r          = 1.1;

k_p_r           = i * omega_n_r^2;
k_d_r           = i * omega_n_r * 2 * zeta_r;

K_p_r           = k_p_r * eye(3);
K_d_r           = k_d_r * eye(3);
K_r             = 1;

torque_control  = @(eta_d, eta, eta_integral)  -K_r * (K_d_r * (eta(4:6)) + K_p_r * (J2' * (eta(1:3) - eta_d(1:3))) );

%Low Pass filtering of Eta_d
T_eta_d         = 0.5;                 %Range (0,0)

%Saturation

sat = @(x, max_ang) min(max(x, -max_ang), max_ang);

% Uncertainties               
if cfg.P_Uncertainty == 1
    q               = 0.03;
else
    q               = 0;
end

Q               = eye(3) * q;

%Systems
x_dot           = @(x,u) (A*x + B*( J1'*u + sqrt(Q) * randn(3,1) ) + [zeros(3,1); m_c * g])/m_c;
omega_dot       = @(omega,M) (I\M);


%% Simulation

alpha           = tau(x,p_integral);
if norm(tau(x,p_integral)) > cfg.k_f / 2
    alpha = alpha/norm(alpha) * (cfg.k_f / 2);
end

u_des           = alpha - m_c * g;

f               = -u_des(3)/(cfg.k_f * cos(eta(2)) * cos(eta(1)));

eta_d_1         = -asin( u_des(2)/(cfg.k_f*f) );
eta_d_2         = asin( u_des(1)/(cfg.k_f*f*cos(eta_d(1))) );

eta_d(1)        = eta_d(1) + T_eta_d*(eta_d_1 - eta_d(1));
eta_d(2)        = eta_d(2) + T_eta_d*(eta_d_2 - eta_d(2));

u               = [ 0; 0; -cfg.k_f * f ];
a               = [zeros(3), eye(3)] * x_dot(x,u);

M               = torque_control(eta_d, eta);

x               = rk4(x_dot, x, u, h);
omega_next      = rk4(omega_dot, eta(4:6), M, h);

eta(1:3)        = eta(1:3) + J2 * eta(4:6) * h;
eta(4:6)        = omega_next;


end