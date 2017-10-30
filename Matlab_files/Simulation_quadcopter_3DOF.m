%% Initialization
helicopter_init();
EKF_init();
%% Simulation

%Variables
x               = zeros(9,1);
x(1:3)          = p_d(:,1);
x_hat           = x ;
P_pri           = 100*eye(9);
eta             = zeros(6,1);
eta_d           = zeros(3,1);
p_integral      = zeros(3,1);
p_d_1_pre       = p_d(:,1);
p_d_pre         = zeros(3,1);

for k = 1:length(time)-1
   
   log.store('p', x(1:3), k);
   log.store('v', x(4:6), k);
   log.store('theta', eta(1:3), k);
   log.store('eta_d', eta_d, k);
   log.store('x_hat', x_hat(1:3), k);
   
   %Low-Pass filter reference signal
   p_d_1_k                  = p_d_1_pre + T_p_d * ( p_d(:,k) - p_d_1_pre);
   p_d_k                    = p_d_pre + T_p_d * ( p_d_1_k - p_d_pre);
   v_d_k                    = calc_diff(p_d_k, p_d_pre, h);
   
   log.store('p_d', p_d_k, k);
   log.store('v_d', v_d_k, k);
   
   %%

   y                                    = h_EK(x, sig_pos) + sqrt(R_PR) * randn(N_measurment,1);
   x_m                                  = x + r_PR*randn(9,1);
   a_m                                  = a + cfg.bias_am;
   
   [P_pri, x_hat, x_prev, error]        = EKF_range(x_hat, sig_pos, y, a_m, P_pri, A_EKF, B_EKF, R_PR, Q_PR, H_EK);
   
   p_integral                           = p_integral + (x_hat(1:3) - p_d(:,k));
   [x_PVA, eta, eta_d, u, M, a]         = Helicopter_3DOF(cfg, x(1:6,1), x_prev(1:6,1), eta, eta_d, p_d_k, v_d_k, p_integral, h);
   x                                    = [x_PVA; x(7:9,1)];
   
   %
   p_d_1_pre                    = p_d_1_k;
   p_d_pre                      = p_d_k;
   
   log.store('error', error, k);
   log.store('x_hat_filtered', x_prev(1:3), k);
   log.store('p_m', x_m(1:3), k);
   log.store('tau', u, k);
   log.store('M', M, k);
   log.store('bias_am', x_hat(7:9), k);
   
   
end

%% Plotting Simulator

p   = log.get('p', ':');
p_d = log.get('p_d', ':');
p_m = log.get('p_m', ':');
figure(1)
plot(time, p);
hold on
plot(time, p_d, '--');
legend('p_x', 'p_y', 'p_z', 'x_d', 'y_d', 'z_d')
hold off
% plot(time, p-p_d);
% legend('e_x', 'e_y', 'e_z')
xlabel ('time \{s\}')
ylabel ('Position \{m\}')
title  ('Position')
grid on

if cfg.Model_3D == 0
    figure(5)
    plot3(p(1,:),p(2,:),p(3,:), 'Linewidth', 4);
    view(45,45);
    cla
    patch([p(1,:) nan],[p(2,:) nan],[p(3,:) nan],[p(3,:) nan],'EdgeColor','interp','FaceColor','none')
    xlabel ('x-axis \{m\}');
    ylabel ('y-axis \{m\}');
    zlabel ('Height \{m\}');
    title  ('Helicopter Path')
    grid on
end

% tau = log.get('tau', ':');
% 
% figure(2)
% plot(time, tau)
% legend('\tau_x', '\tau_y', '\tau_z')
% xlabel ('time')
% ylabel ('Controller gain')
% title  ('Controller')
% grid on
% 
% M = log.get('M', ':');

% figure(4)
% plot(time, M)
% legend('M_\phi', 'M_\theta', 'M_\psi')
% xlabel ('time')
% ylabel ('Angular Controller gain')
% title  ('Angular Controller')
% grid on

theta = log.get('theta', ':');
eta_d = log.get('eta_d', ':');

figure(3)
plot(time, theta)
hold on
plot(time, eta_d, '--');
hold off
legend('\phi', '\theta', '\psi', '\phi_d', '\theta_d', '\psi_d')
xlabel ('time')
ylabel ('Angles')
title  ('Attitude')
grid on

%% Plotting EKF

%Position Estimate
p       = log.get('p', ':');
p_hat_m = log.get('x_hat', ':');
p_hat_p = log.get('x_hat_filtered', ':');

figure(5)
plot(time, p(1,:));
hold on
plot(time, p_hat_m(1,:), '--');
plot(time, p_hat_p(1,:), '--');
%plot(time, p_m(1,:), '-.');
legend({'$p_x$', '$\hat{x}^-$', '$\hat{x}^+$'}, 'Interpreter', 'latex')
hold off
xlabel ('time \{s\}')
ylabel ('Position \{m\}')
title  ('Position')
grid on

%Bias estimation
bias_am = log.get('bias_am', ':');
figure(6)
plot(time, low_pass_filter(bias_am,400));
hold on
plot(time, ones(3,length(bias_am)).*cfg.bias_am);
legend({'$\hat{\beta}_x$', '$\hat{\beta}_y$', '$\hat{\beta}_z$', '$\beta_x$', '$\beta_y$' '$\beta_z$'}, 'Interpreter', 'latex', 'FontSize',16)

error = log.get('error', ':');
checkNormalDistribution(7,error(:), '$\varepsilon_y$');

%% Save Log
if cfg.SaveLog == 1
    storeName = '_Simulation';
    TimeNow = datestr(now, 'yyyy.mm.dd-HH.MM.SS');
    save(['logs/' 'sim_' storeName TimeNow '.mat'], 'log');
end
