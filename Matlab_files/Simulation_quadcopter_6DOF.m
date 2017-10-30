%% Initialization
functions;

%% Switches

%% Parameters

% Simulation
h               = 0.05;
tend            = 20;

time            = 0:h:tend;

% Regulator
K_d             = 1*diag(3);
K_p             = 0.5*diag(3);

k_f             = 1;

max_ang         = pi/4;

% Physical
m_c             = 5;
g_z             = 9.81;

g               = [0; 0; g_z];
w               = zeros(3,1);

I               = [1, 0, 0;
                   0, 1, 0;
                   0, 0, 1];
          

%% Variables
% References

log = Logger();
log.init(time,h);

%DOF variables and their derivates
log.add('p',3);
log.add('v',3);
log.add('theta',3);
log.add('omega',3);

%Regulator Variables
log.add('tau',3);
log.add('M',3);

%References
log.add('p_d',3);
log.add('v_d',3);
log.add('theta_d',3);


%% Simulation

ref             = @(T) [cos(T); sin(T); 1];

p_d             = [cos(t/4)+1; sin(t/4)+1; ones(size(t))];
v_d             = diff(p_d,1,2)./h;

log.set('p_d', p_d);
log.set('v_d', [v_d, zeros(3,1)]);

bound_desired_angle     = @(theta) mod(theta,pi) - pi/2;
sat                     = @(x,min_x,max_x) min(max_x,max(min_x,x));

for k = 2:length(time)
    
    
    
    t = time(k);
    tau_d               =  -m_c * K_d * ( v(:,t) - v_d(:,t) ) - m_c * K_p * ( p(:,t) - p_d(:,t) );
    f                   = ( -m_c * g(3) + tau_d(3,t) ) / (k_f * cos( theta(2,t) ) * cos( theta(1,t) ) );
    
    theta_d(1,t)        = bound_desired_angle(asin(sat( tau_d(2,t)/(k_f * f),                  -sin(max_ang),sin(max_ang))));
    theta_d(2,t)        = bound_desired_angle(asin(sat(-tau_d(1,t)/(k_f*f * cos(theta_d(1,t))),-sin(max_ang),sin(max_ang))));
    
    [~,~,T]             = eulerang(theta(1,t), theta(2,t), theta(3,t));
    M(:,t)              = -K_d * omega(:,t) - T' * K_p * (theta(:,t) - theta_d(:,t));
    
    tau(:,t)            = (-k_f*f).*[ sin(theta(3,t))*sin(theta(1,t)) + cos(theta(3,t))*cos(theta(1,t))*sin(theta(2,t));
                                     -cos(theta(3,t))*sin(theta(1,t)) + sin(theta(2,t))*sin(theta(3,t))*cos(theta(1,t));
                                                           cos(theta(2,t))*cos(theta(1,t))                             ];
    
    v_dot(:,t)          = g + R_theta(theta(1,t), theta(2,t), theta(3,t))*tau(:,t)/m_c + w/m_c;
    theta_dot(:,t)      = T*omega(:,t);
    omega_dot(:,t)      = I\S(I*omega(:,t))*omega(:,t) + I\M(:,t);
    
    
    %Updates
    p(:,t+1)            = p(:,t)     + h*v(:,t);
    v(:,t+1)            = v(:,t)     + h*v_dot(:,t);
    omega(:,t+1)        = omega(:,t) + h*omega_dot(:,t);
    theta(:,t+1)        = theta(:,t) + h*theta_dot(:,t);
end

%% Plotting

figure
plot3(p(1,:),p(2,:),p(3,:), 'Linewidth', 4);
cla
patch([p(1,:) nan],[p(2,:) nan],[p(3,:) nan],[p(3,:) nan],'EdgeColor','interp','FaceColor','none')
xlabel ('x-axis');
ylabel ('y-axis');
zlabel ('Height');
title  ('Helicopter Path')
grid on