%% This file simulates and tries to control the suspended pendulum system in 5DOF

clear; 

addpath('helpers');
% Set up some parameters
pL = 4;
pm_c = 2.3;
pm_L = 0.450;
pg = 9.81;
pd = 0.01;
epsilon = 0.001;

% Set Case Number!
cTracking = Case.get('Tracking');
cBox      = Case.get('Box');
cBoxShape = Case.get('BoxShape');
cLineFeedback = Case.get('LineFeedback');

cStep     = Case.get('Step');
cZigZag   = Case.get('ZigZag');

csNone = SubCase('slung');
csInput = SubCase('slung,input');
csDelayed = SubCase('slung,delayed');
csBoth    = SubCase('slung,input,delayed');
csBothWind    = SubCase('slung,input,delayed,wind');
csPid   = SubCase('');


% which to run?
doCase = cZigZag;
doSubCase = SubCase('slung,wind,input,delayed');

% Just a backup
if ~isa(doSubCase, 'SubCase')
    error('Wrong sub case specified. Using default. ');
    doSubCase = SubCase();
end

% Activate the plot at all
doInteractivePlot = 1;
doOnlyLastFrame = 1;


% Simulation Options
tend = 10;
h = 0.05;

% Disturbance
bias = 0*[1.4 0.5 0.1]';


wind = 0*[3 3 0]';

pm_L_shaper = pm_L*1;
pL_shaper = pL*1;


% Control parameters
k1 = 3.397;
k2 = 12.6583;
gamma_m = 5;
K_i = 3;

k1 = 0.4226;
k2 = 4.7321;
K_i = 1;

k2_45 = 4*k2;

% Try to keep actual Ki the same, which is
% Ki = rho * k1;

rho = 1.5774 / 1;

K_i = 1.5774;

% Try two, same dampings etc
%k1 = 1.5774;
%k2 = 1.2679;


k2_45 = 4*k2;

% Reference model parameters
refmodel_max_v = 4;
refmodel_max_a = 7;
refmodel_max_j = 6;
refmodel_w0 = 1.1;
refmodel_zeta = 0.9;

% Control program flow
constantRef = 1;

startAltitude = 0;
endAltitude = -2;
timeToMaxAltitude = tend-8;

% Use adaptive controller
useAdaptiveController = 0;
testAdaptiveM = 1;

% Check if use reference model
useReferenceModel = 1;

% Use input shaper?
useInputShaperZVD = 0;
useInputShaperZV  = 1;

% Do Wps. % Only used in constant ref!
doWP = 1;
WPs = 10*[0 1 0;
       1 1 -1; 
       1 0 -1; 
       0 0 0;]';

WPsWaitFor = [];
exactAlphaIntegration = 0;
useOde45              = 0;

% Set negative to not activate. 
accelRefFeedbackThreshold = -0.05;
   
   
if doCase == cTracking
    constantRef = 0;
    useReferenceModel = 0;
    % Use input shaper?
    useInputShaperZVD = 0;
    useInputShaperZV  = 0;
    
    tend = 10;
    timeToMaxAltitude = tend-8;
    timeActivateFeedback = -1;
    
elseif doCase < cStep
    constantRef = 1;
    useReferenceModel = 1;
    
    tend = 40;
    timeActivateFeedback =0.001;
    
    % Use input shaper?
    useInputShaperZV  = 1;
    if doCase == cBoxShape
        useInputShaperZVD = 0;
        useInputShaperZV  = 1;
    end
    
    if doCase == cLineFeedback
        useReferenceModel = 1;
        doWP = 0;
        tend = 16;
        useInputShaperZV  = 0;
        timeActivateFeedback = 9.5;
        timeActivateFeedback = 5.5;
        timeActivateFeedback = 0;
    end
end

if doCase == cStep || doCase == cZigZag
    constantRef = 1;
    useReferenceModel = 1;
    doWP = 1;
    useInputShaperZV  = 0;
end

if doCase == cStep
    tend = 70;
    WPs = Rzyx(0,0,deg2rad(67)) *[ 0  0 0;
            20 0 0;
            0  0 0]';
    
    WPsWaitFor = [18 30];   
elseif doCase == cZigZag
    tend = 80;
    WPs = Rzyx(0,0,deg2rad(67)) *[0  0  0;
               10 0  0;
               20 0  0; 
               0  20 0; 
               20 20 0;
               0  0  0;
               10 0  0]';
    WPsWaitFor = [0 37.5];       
           
end

wind = [0 0 0]';
if doSubCase.enWind 
    wind = [1 -3 0]';
end

if doSubCase.enInput
   useInputShaperZV  = 1;
end

accelRefFeedbackThreshold = -0.05;
if doSubCase.enDelayed
   timeActivateFeedback =0.001;
   accelRefFeedbackThreshold = 0.05;
end
    

% Initialize anon functions
addpath('generated/');
M = @(eta)(f_5dof_MassMatrix(eta, pL, pm_c, pm_L));
C = @(eta, nu)(f_5dof_CoreolisMatrix(eta, nu, pL, pm_c, pm_L));
G = @(eta)(f_5dof_gravity(eta, pL, pm_c, pm_L, pg));
D = diag([0 0 0 pd pd]);
H = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0];
 
D_13 = diag([0.15 0.15 0.2 0 0]);
 
% Restate for singularity avoidance
M = @(eta)(f_5dof_MassMatrix_singularity_avoidance(eta, pL, pm_c, pm_L, epsilon));
C = @(eta, nu)(f_5dof_CoreolisMatrix_singularity_avoidance(eta, nu, pL, pm_c, pm_L, epsilon));
 
 
M_hat = @(eta, m_L)(f_5dof_MassMatrix(eta, pL, pm_c, m_L));
C_hat = @(eta, nu, m_L)(f_5dof_CoreolisMatrix(eta, nu, pL, pm_c, m_L));
G_hat = @(eta, m_L)(f_5dof_gravity(eta, pL, pm_c, m_L, pg));
 
M_L = @(eta)(f_5dof_MassMatrix(eta, pL, 0, 1));
C_L = @(eta, nu)(f_5dof_CoreolisMatrix(eta, nu, pL, 0, 1));
G_L = @(eta)(f_5dof_gravity(eta, pL, 0, 1, pg));


 
K1 = diag([k1 k1 k1]);
K2 = diag([k2 k2 k2 k2_45 k2_45]);
 
% Simulation options



time = 0:h:tend;

% Initial conditions
eta_0 = [0, 0, 0, 0.00, 0.00]';

if doCase == Case.get('zigzag')
    eta_0(1:3) = WPs(:,2);
end

nu_0  = [0, 0, 0, 0, 0]';

% Desired positions.

% y = 1.5;
desiredPos = [-4,-4, 0]';

r   = zeros(3, length(time));
dr  = zeros(3, length(time));
d2r = zeros(3, length(time));
d3r = zeros(3, length(time));

r(:,1) = eta_0(1:3);


if constantRef && ~doWP
   r = repmat(desiredPos, 1, length(time));
   
   
   % Try to do convolution!
   

   
   
elseif ~constantRef
   f = 1.5;
   
   f = 1;
   
   altDiff = endAltitude - startAltitude;
   altFreq = altDiff/timeToMaxAltitude;
   
   % Create a nice spiral
   % Get altitude theta, as time up untin timeToMaxAltitude
   % Get Index:
   altMaxTimeIndex = find(time == timeToMaxAltitude);
   altTime = [time(1:altMaxTimeIndex) time(altMaxTimeIndex)*ones(1, length(time)-altMaxTimeIndex)];
   
   %r =  2*pL*[cos(f*time); sin(f*time); altFreq*altTime];
   r =  2*pL*[cos(f*time); sin(f*time); 0*-0.1*f*time];
   
   
   % Change so that we go to senter pos after that time
   %r = [r(1:2, 1:altMaxTimeIndex) zeros(2, length(time)-altMaxTimeIndex);
   %     r(3, :)];
   
   r(1,:) = eta_0(1) + r(1,:) - r(1,1);
   r(2,:) = eta_0(2) + r(2,:) - r(2,1);
   r(3,:) = eta_0(3) + r(3,:) - r(3,1);
   
   dr = 2*pL*f*[-sin(f*time); cos(f*time); zeros(1, length(time))];
   d2r = 2*pL*f^2*[-cos(f*time); -sin(f*time); zeros(1, length(time))];
   d3r = 2*pL*f^3*[sin(f*time); -cos(f*time); zeros(1, length(time))];
end



% Initialize logger

log = Logger();

log.init(time, h);
log.add('Eta', 5);
log.add('Nu', 5);
log.add('tau', 5);
log.add('tau_45', 2);
log.add('z2_integral', 3);
log.add('dNu', 5);
log.add('alpha', 5);
log.add('m_L_hat', 1);
log.add('r', 3);
log.add('dr', 3);
log.add('d2r', 3);
log.add('r_preconv',3);
log.add('dr_preconv', 3);
log.add('d2r_preconv', 3);
log.add('r_prerefmodel', 3);
log.add('alpha_45', 2);
log.add('tau_g', 3);
log.add('tau_beta_1', 3);
log.add('tau_beta_2', 3);
log.add('tau_main', 3);
log.add('feedback_active', 1);
log.add('r_postconv', 3);
log.add('dr_postconv', 3);
log.add('d2r_postconv', 3);

% Store that reference before ref. model
log.set('r_prerefmodel', r);
% some options
m_L_min = 0.001;
m_L_max = 10;

% start simulation loop
eta = eta_0;
nu  = nu_0;
alpha_45 = [0; 0];
m_L_hat = 0.001;

% Store initial position
log.store('Eta', eta, 1);


% Use reference model in stead of reference derivatives

refModelState = [eta_0(1:3) ; zeros(3,1); zeros(3,1); zeros(3,1)];

refModelState_x = [eta_0(1); 0; 0];
refModelState_y = [eta_0(2); 0; 0];
refModelState_z = [eta_0(3); 0; 0];

fprintf('*** Starting Simulation ***\n');
c_sec = tic;
c_total = tic;
reverseString = '';
step = 0;
t_sec = 0;


       
curWP = 1;
near = 0;
near_time = 0;
for k = 2:length(time)
   t = time(k);
   
   
   
   
   % Check if we should apply the reference model
   if useReferenceModel
        w0=2; % this a go starting point increase for faster decrease for slower
        zeta=1/sqrt(2);
        %zeta=0.9;
        
        I = eye(3);
        z = zeros(3);

        Add=[   z I z ;
                z z I;
                -w0^3*I -(2*zeta+1)*w0^2*I -(2*zeta+1)*w0*I];
        Bdd=[z; z; w0^3*I];
        Cdd=eye(9);
        Ddd=[ z z z]';
        
        f = @(x, u)(Add*x + Bdd*u);
        
        
        % Do reference model properly!
        
        % Set parameters
        w0 = refmodel_w0;
        xi = refmodel_zeta;
        
        % State rename
        x = refModelState;
        p = x(1:3);
        v = x(4:6);
        a = x(7:9);
        j = x(10:12);
        
        % Set reference
        if doWP && constantRef
            ref = WPs(:, curWP);
        else
            ref = r(:,k);
        end
        
        k3_ref = (2*xi+1)*w0;
        k2_ref = k3_ref\(2*xi + 1)*w0^2;
        k1_ref = (k2_ref*k3_ref)\w0^3;
        
        k4_ref = (4*xi)*w0;
        k3_ref = k4_ref\(2 + 4*xi^2)*w0^2;
        k2_ref = (k4_ref*k3_ref)\((4*xi)*w0^3);
        k1_ref = (k4_ref*k3_ref*k2_ref)\w0^4;

        % Step 1: v-controller
        tau1 = k1_ref * (ref - p);

        if norm(tau1) > refmodel_max_v
            tau1 = refmodel_max_v * tau1 / norm(tau1);
        end

        % Step 2: a-controller
        tau2 = k2_ref * (tau1 - v);

        if norm(tau2) > refmodel_max_a
            tau2 = refmodel_max_a * tau2 / norm(tau2);
        end

        % Step 3: j-controller
        tau3 = k3_ref * (tau2 - a);
        if norm(tau3) > refmodel_max_j
            tau3 = refmodel_max_j * tau3 / norm(tau3);
        end

        % Step 4: dj controller
        tau4 = k4_ref * (tau3 - j);

        % Integrate
        x = x + h*[x(4:12); tau4];

        % Restate
        refModelState = x;

        
        % Integrate 
        
        if doWP && constantRef
            %refModelState_x = rk4(f, refModelState_x, WPs(1, curWP), h);
            %refModelState_y = rk4(f, refModelState_y, WPs(2, curWP), h);
            %refModelState_z = rk4(f, refModelState_z, WPs(3, curWP), h);
            
            %refModelState = rk4(f, refModelState, WPs(:, curWP), h);
        else
            %refModelState_x = rk4(f, refModelState_x, r(1, k), h);
            %refModelState_y = rk4(f, refModelState_y, r(2, k), h);
            %refModelState_z = rk4(f, refModelState_z, r(3, k), h); 
            
            %refModelState = rk4(f, refModelState, r(:,k), h);
        end
        
        if norm(refModelState(4:6)) > 1.1*refmodel_max_v
            %refModelState(4:6) = 3*refModelState(4:6)/norm(refModelState(4:6));
            warning('Should not happen. V max');
        end
        
        if norm(refModelState(7:9)) > 1.1*refmodel_max_a
            %refModelState(7:9) = 3*refModelState(7:9)/norm(refModelState(7:9));
            warning('Should not happen. a max');
        end
        
        % Update desired positions
        %%2,k) = refModelState_y(1);
        %r(3,k) = refModelState_z(1);
        r(:,k) = refModelState(1:3);
        
        %dr(1,k) = refModelState_x(2);
        %dr(2,k) = refModelState_y(2);
        %dr(3,k) = refModelState_z(2);
        dr(:,k) = refModelState(4:6);

        
        %d2r(1,k) = refModelState_x(3);
        %d2r(2,k) = refModelState_y(3);
        %d2r(3,k) = refModelState_z(3);
        d2r(:,k) = refModelState(7:9);
        
        d3r(:,k) = refModelState(10:12);

   end
   
   if norm(r(:,k) - WPs(:,curWP)) < 1.5*refmodel_max_v % && curWP < size(WPs, 2)
       if ~near 
           near = 1;
           near_time = time(k);
       end
       
       if near && (length(WPsWaitFor) < curWP || time(k) - near_time > WPsWaitFor(curWP))
         curWP = curWP + 1;
         near = 0;
       end
       if curWP > size(WPs, 2)
           curWP = size(WPs, 2);
       end
       k;
   end
   
end
% Log 
log.set('r_preconv', r);
log.set('dr_preconv', dr);
log.set('d2r_preconv', d2r);

if useInputShaperZV || useInputShaperZVD
    % Create shaper
    omega_n = sqrt(pg/pL_shaper);

    xi = 0.2; % this is just a guess!
    xi = pd/(2*omega_n*pm_L_shaper);

    omega_d = omega_n*sqrt(1-xi^2);
    

    Td = 2*pi/omega_d;

    K = exp(-(xi*pi)/sqrt(1-xi^2));
end

% Apply shaper
if useInputShaperZVD
    

    A1 = 1/(1 + 2*K + K^2);
    A2 = 2*K/(1 + 2*K + K^2);
    A3 = K^2/(1 + 2*K + K^2);

    t2 = Td/2;
    t3 = Td;

    % Then, create a sequence at these times. 

    % signal length:
    k3 = find(time >= t3,1);
    k2 = find(time >= t2,1);

    shaper = zeros(1, k3);
    shaper(k3) = A3;
    shaper(k2) = A2;
    shaper(1)  = A1;   
    
    
    % Try to "up" the shaper by initial condition
    r_n = [];
    r_n(1,:) = conv(r(1,:), shaper+ r(1,1)*ones(size(shaper)));
    r_n(2,:) = conv(r(2,:), shaper+ r(2,1)*ones(size(shaper)));
    r_n(3,:) = conv(r(3,:), shaper+ r(3,1)*ones(size(shaper)));
    %r = r_n;
    
    % Does not work as expected, using zero initialized version
    r = convn(r, shaper);
    dr = convn(dr, shaper);
    d2r = convn(d2r, shaper);
    
elseif useInputShaperZV
    
    A1 = 1/(1 + K);
    A2 = K/(1 + K);

    t2 = Td/2;
    

    % Then, create a sequence at these times. 

    % signal length:
   
    k2 = find(time >= t2,1);

    shaper = zeros(1, k2);
    shaper(k2) = A2;
    shaper(1)  = A1;   
    
    
    % Try to "up" the shaper by initial condition
    r_n = [];
    r_n(1,:) = conv(r(1,:), shaper+ r(1,1)*ones(size(shaper)));
    r_n(2,:) = conv(r(2,:), shaper+ r(2,1)*ones(size(shaper)));
    r_n(3,:) = conv(r(3,:), shaper+ r(3,1)*ones(size(shaper)));
    %r = r_n;
    r = convn(r, shaper);
    dr = convn(dr, shaper);
    d2r = convn(d2r, shaper);
    
    % Add initial position to the first k < k2 
    r(:,1:k2-1) = r(:,1:k2-1) + A2*eta_0(1:3)*ones(1,k2-1);
end

log.set('r', r(:,1:length(time)));
log.set('r_postconv', r(:,1:length(time)));
log.set('dr_postconv', dr(:,1:length(time)));
log.set('d2r_postconv', d2r(:,1:length(time)));


doAnalysis = 0;

if doAnalysis
    profile on
end



% Start dynamic simulation
z2_integral = zeros(3, 1);
%z2_integral = wind.*diag(D_13(1:3,1:3));
lookup = @(A, index)(subsref(A, struct('type', '()', 'subs', {index})));
for k = 1:length(time)-1
   t = time(k);
   refK = k;
   
   if doAnalysis && k > 2/h
       break
   end
   
   % Do some control
   tau = zeros(5,1); 
   tau(3) = - pg*(pm_L + pm_c);
   

   % Ait, check feedback stuff
   nSamplesBackwards = floor(6/h);
   sigmoid = @(t)(1./(1+exp(-t)));
   
   activeCouseCase = (doCase == cLineFeedback && t > timeActivateFeedback);
   activeCouseTime = (doCase == cBox || doCase == cBoxShape) && ...
                        (timeActivateFeedback > 0 && t > timeActivateFeedback);
   
   activeCouseAccel = refK > nSamplesBackwards ...
                       && accelRefFeedbackThreshold > 0;
                   
     %all(colnorm(log.get('d2r_postconv',:,refK-nSamplesBackwards:refK)) < accelRefFeedbackThreshold ))
   activeGain = 0;
   
   
   if activeCouseCase || activeCouseTime
       activeGain = 1;
   end
                    
   if activeCouseCase || activeCouseTime || activeCouseAccel
            
        omega_n = sqrt(pg/pL);

        xi = 0.2; % this is just a guess!
        xi = pd/(2*omega_n*pm_L);
        omega_d = omega_n*sqrt(1-xi^2);
        Td = 2*pi/omega_d;
        tau_d = 0.325*Td;
        Gd    = 0.325;
        oldAngles = [0; 0];
        olddAngles = [0; 0];
        oldd2Angles = [0; 0];
        
        % Gain activation test
        if activeCouseAccel
            % Percentage og samples below threshold
            pbt = (nSamplesBackwards+1)\sum(colnorm(log.get('d2r_postconv',':',refK-nSamplesBackwards:refK)) < accelRefFeedbackThreshold );
            
            % a = 0:0.01:1;
            % figure(2); clf; plot(a, sigmoid(20*(-0.7+a))); grid on
            activeGain = sigmoid(30*(-0.7+pbt));
            %activeGain = sigmoid(15*(-0.5+pbt));
        else
            activeGain = 1;
        end
        
        Gd = activeGain*Gd;
        
        if t > tau_d
            % Find time
            oldTimeIndex = find(time>t-tau_d,1,'first');
            oldAngles = log.get('Eta', 4:5, oldTimeIndex);
            olddAngles = log.get('Nu', 4:5, oldTimeIndex);
            oldd2Angles = log.get('dNu', 4:5, oldTimeIndex);
        end
        
        %Re-define current desired pos. 
        log.store('feedback_active', activeGain, k);
        k_ = refK;
        r(1,k_) = r(1,k_) + Gd*pL*sin(oldAngles(2));
        r(2,k_) = r(2,k_) - Gd*pL*sin(oldAngles(1));
        
        dr(1,k_) = dr(1,k_) + Gd*pL*cos(oldAngles(2))*olddAngles(2);
        dr(2,k_) = dr(2,k_) - Gd*pL*cos(oldAngles(1))*olddAngles(1);
        
        d2r(1,k_) = d2r(1,k_) + Gd*pL*(-sin(oldAngles(2))*olddAngles(2) + cos(oldAngles(2))*oldd2Angles(2));
        d2r(2,k_) = d2r(2,k_) - Gd*pL*(-sin(oldAngles(1))*olddAngles(1) + cos(oldAngles(1))*oldd2Angles(1));
   end
   
   % Define error variables
   z1 = H*eta - r(:,refK);
   
   % Virtual control 1-2
   alpha_13 = dr(:,refK) - K1 * z1;
   alpha = [alpha_13; alpha_45];
   
   z2 = nu - alpha;
   
   dalpha_13 = d2r(:,refK) - K1*H*z2 + K1*K1*z1;
   dalpha_13_f = @(eta, nu)(d2r(:,refK) - K1*(nu(1:3)-dr(:,refK)));
   
   % Shorthands
   z45 = z2(4:5);
   theta_L = eta(5);
   phi_L = eta(4);
   
   G45 = G(eta);
   G45 = G45(4:5);
   G45_f = @(eta)(lookup(G(eta), {4:5}));
   
   K45 = K2(4:5,4:5);
   
   M4513 = M(eta);
   M4513 = M4513(4:5,1:3);
   M4513_f = @(eta)(lookup(M(eta), {4:5, 1:3}));
   M1345_f = @(eta)(lookup(M(eta), {1:3, 4:5}));
   
   Malpha = M(eta);
   Malpha = Malpha(4:5, 4:5);
   Malpha_f = @(eta)(lookup(M(eta), {4:5, 4:5}));
   
   Calpha = C(eta, nu);
   Calpha = Calpha(4:5, 4:5);
   Calpha_f = @(eta, nu)(lookup(C(eta, nu), {4:5, 4:5}));
   
   C1345_f = @(eta, nu)(lookup(C(eta, nu), {1:3,4:5}));
   Dalpha = D(4:5, 4:5);
   
   % Define residual
   gamma = @(eta, pm_L)(-G45 + K45*z45 - M4513*dalpha_13);
   gamma_f = @(eta, nu, alpha_45)(-G45_f(eta) + K45*(nu(4:5)-alpha_45) - M4513_f(eta)*dalpha_13_f(eta,nu));
   
   
   
   dalpha_45 = @(alpha_45, ~)( Malpha\(-Dalpha*alpha_45 - Calpha*alpha_45 + gamma(alpha_45, pm_L)));
   dalpha_45_f = @(eta, nu, alpha_45)( Malpha_f(eta)\(-Dalpha*alpha_45 - Calpha_f(eta, nu)*alpha_45 + gamma_f(eta, nu, alpha_45)));
   
   
   dalpha = [dalpha_13; dalpha_45(alpha_45, 0)];
   
   
   % Integrate z2
   %z2_integral = z2_integral + h * z1;

   z2_integral = z2_integral + h * K_i * H*z2;
   %z2_integral = z2_integral + h * K_i * z1;
   %z2_integral = 0*z2_integral;
   % Main control part
   tau = C(eta, nu)*alpha + D*alpha + G(eta) - H'*z1 + M(eta)*dalpha - K2*z2 - H'*z2_integral + D_13*nu;
   tau_f = @(eta, nu, alpha_45, t)(H'*C1345_f(eta, nu)*alpha_45 + H'*M1345_f(eta)*dalpha_45_f(eta, nu, alpha_45) + G(eta) + ...
       +M(eta)*H'*d2r(:,refK) - H'*(1+K2(1:3,1:3)*K1)*(eta(1:3)-r(:,refK)) -H'*(H*M(eta)*H'*K1+K2(1:3,1:3))*H*H'*(nu(1:3)-dr(:,refK)));
   
   % Redo, for some logging
   M_ = M(eta);
   C_ = C(eta, nu);
   
   tau_g = G(eta);
   tau_beta_1 = C_(1:5,4:5)*alpha(4:5,1) + D*alpha;
   tau_beta_2 = M_(1:5,4:5)*dalpha(4:5,1);
   tau_main = -H'*z1 + M_(1:5,1:3)*dalpha(1:3,1) -K2*z2 - H'*z2_integral + D_13*nu;
   % For fun, print third value
   %tau(3)
   
   if ~doSubCase.enSlung
       tau = tau_g + tau_main;
   end
   
   log.store('tau_g', tau_g(1:3),k);
   log.store('tau_beta_1', tau_beta_1(1:3), k);
   log.store('tau_beta_2', tau_beta_2(1:3), k);
   log.store('tau_main', tau_main(1:3), k);
   log.store('z2_integral', z2_integral, k);
   

   if useAdaptiveController || testAdaptiveM
       
       % Main control part
       if useAdaptiveController
        tau = C_hat(eta, nu, m_L_hat)*alpha + G_hat(eta, m_L_hat) - H'*z1 + M_hat(eta, m_L_hat)*dalpha - K2*z2 + D*nu - H'*z2_integral;
       end
      
       % Integrate estimate
       Phi = M_L(eta)*dalpha + C_L(eta, nu)*alpha + G_L(eta);
       dm_L_hat = @(~, ~)(- gamma_m * z2' * Phi);
       m_L_hat_next = rk4(dm_L_hat, 0, 0, h);
       
       % Hackish-projection function
       if m_L_hat_next < m_L_min
           m_L_hat_next = m_L_min;
       elseif m_L_hat_next > m_L_max
           m_L_hat_next = m_L_max;
       end
       
       m_L_hat = m_L_hat_next;
       
       % Log
       log.store('m_L_hat', m_L_hat, k);
       
   end
   
   % Log third value
   log.store('tau_45', tau(4:5), k);
   
   % aand, set to zero
   tau(4:5) = [0; 0];
   
   
   if t > tend/2
       %tau(1:2) = [0; 0];
   end
   
   
   
   
   f_eta = @(eta, nu, tau)(nu);
   f_nu = @(eta, nu, tau)(M(eta)\(H'*bias + tau - C(eta, nu)*nu - G(eta) - D*nu - D_13*(nu - H'*wind)));
   f_alpha_45 = @(eta, nu, alpha_45)(dalpha_45_f(eta, nu, alpha_45));
   dnu = f_nu(eta, nu, tau);
   

   
   if exactAlphaIntegration
       f = @(x, tau,t)([f_eta(      x(1:5), x(6:10), tau); ...
                      f_nu(      x(1:5), x(6:10), tau);
                      f_alpha_45(x(1:5), x(6:10), x(11:12))]);
       f_tau = @(x, t)(H'*H*tau_f(x(1:5), x(6:10), x(11:12), t));
       
       if useOde45
           f_ode = @(t, x)(f(x, f_tau(x, t),t));
           %nu_next  = rk4(f, nu, tau, h);
           
           options = odeset('RelTol',1e-4);

           res = ode15s(f_ode, [0 h], [eta; nu; alpha_45], options);
           x_next = res.y(:,end);
       else
           x_next = rk4(f, [eta; nu; alpha_45], f_tau, h, 0);
       end
       eta = x_next(1:5);
       nu  = x_next(6:10);
       alpha_45 = x_next(11:12);
   else
       f = @(x, tau)([f_eta(x(1:5), x(6:10), tau); ...
                      f_nu( x(1:5), x(6:10), tau)]);
       
       
       %nu_next  = rk4(f, nu, tau, h);
       x_next = rk4(f, [eta; nu], tau, h);
       
       % Integrate to get alpha_3
        alpha_45 = rk4(dalpha_45, alpha_45, 0, h);

        % Update current state
        %eta = eta_next;
        %nu  = nu_next;

        eta = x_next(1:5);
        nu  = x_next(6:10);
   end
   

   
  
   
   %f = @(nu, tau)(M(eta)\(H'*b + tau - C(eta, nu)*nu - G(eta) - D*diag(abs(nu))*nu));

   
   % Combine the integration steps, do it properly
   
   
   
   % Log
   log.store('Eta', eta, k+1);
   log.store('Nu', nu, k+1);
   log.store('dNu', dnu, k+1);
   log.store('tau', tau, k);
   log.store('alpha_45', alpha_45, k+1);

   
   % Print
   checkEveryInterval = 0.4;
    if ~mod(time(k), checkEveryInterval) 
        %time(t)
        step = time(k);
        t_sec = toc(c_sec)/checkEveryInterval;
        if time(k) == 0
            t_sec = 1;
        end
        %t_sec = toc(c_total);
        
        c_sec = tic;
    end
    
    outline = sprintf(['* Simulation time: %.2fs of %ds. \n' ...
                       '* Computational time: %.2fs. \n' ...
                       '* Time remaining: %.1fs. \n'], time(k), tend, toc(c_total), abs(t_sec*tend - toc(c_total))); 
    fprintf([reverseString, outline]);
    reverseString = repmat(sprintf('\b'), 1, length(outline));

   
   
end

% Re-set r log (in case of delayed)
log.set('r', r(:,1:length(time)));
log.set('dr', dr(:,1:length(time)));
log.set('d2r', d2r(:,1:length(time)));

% double-store the last integral value
log.store('z2_integral', log.get('z2_integral', ':', length(time)-1), length(time));
%% Print Done!
t_total = toc(c_total);
outline = sprintf(['*** Simulation complete *** \n' ...
         '* Total computation time: %f \n' ...
         '* Average time per second: %f \n' ...
         '* System Equations: %d \n'...
         '***\n' ...
         ], t_total, t_total/tend, 5);
fprintf([reverseString, outline]);
reverseString = repmat(sprintf('\b'), 1, length(outline));

if doAnalysis
    profile viewer
    return
end

%% Store logs

switch doCase
    case cTracking
        log_cTracking = log;
    case cBox
        log_cBox = log;
    case cBoxShape
        log_cBoxShape = log;
    case cLineFeedback
        log_cLineFeedback = log;
    case cStep
        log_cStep = log;
end


%% Save logs
timeNow = datestr(now, 'yyyy.mm.dd-HH.MM.SS');

% save the log!
switch doCase
    case cTracking
        storeName = 'Tracking_';
        logName = 'log_cTracking';
    case cBox
        storeName = 'Box_';
        logName = 'log_cBox';
    case cBoxShape
        storeName = 'BoxShape_';
        logName = 'log_cBoxShape';
    case cLineFeedback
        storeName = 'LineFeedback_';
        logName = 'log_cLineFeedback';
    case {cStep, cZigZag}
        storeName = [Case.get(doCase) '_' doSubCase.getString() '_'];
        logName = 'log';
end

save(['logs/' 'sim_' storeName timeNow '.mat'], logName);

%% PLOT
% Do some simple plotting!
% figure(2); clf
% 
% subplot(3,1,1);
% title('x');
% plot(time, log.get('Eta', 1));
% 
% subplot(3,1,2);
% title('y');
% plot(time, log.get('Eta', 2));
% 
% subplot(3,1,3);
% title('\theta_L');
% plot(time, log.get('Eta', 3));



% Activate super-screen? 
activateSuperScreen = 1;





p = Plotter;
p.setSize(2,3);
p.activate('main');
%p.activate('theta');
p.activate('pos_error');
%p.activate('alpha_danglular');
p.activate('loadEnergy');
p.activate('z2_integral');


p.activate('vel');

p.activate('vel_error');
p.activate('acc_error');

p.activate('tau');
p.activate('pos');
p.activate('m_L_hat');
% p.activate('loadConstraint');
% p.activate('mu');





% Prepare some data
load = 1;

% Calculate energy of suspended load
loadEnergy = zeros(2, length(time));

for k = 1:length(time)
    
    eta = log.get('Eta', ':', k);
    nu  = log.get('Nu', ':', k);
    
   loadEnergy(1,k) = nu'*M(eta)*nu + pm_L*pg*pL*(1-cos(eta(4))*cos(eta(5)));
   
   % Try to extract only pendulum
   M_now = M(eta);
   eta_ang = eta(4:5);
   nu_ang  = nu(4:5);
   
   loadEnergy(2,k) = nu_ang'*M_now(4:5, 4:5)*nu_ang + pm_L*pg*pL*(1-cos(eta(4))*cos(eta(5)));
   
end



if doInteractivePlot


fig = figure(1); clf; %subplot(3,3,[1 2 4 5 7 8]); hold on;



[on, subcell] = p.isActive('main');
if on
    subplot(subcell{1}{:})
    view(3);
    set(gca, 'ZDir', 'reverse');
    set(gca, 'YDir', 'reverse');
    hold on;
    
    copter_x = log.get('Eta', 1);
    copter_y = log.get('Eta', 2);
    copter_z = log.get('Eta', 3);
    lim_x = [min(copter_x)-1.1*pL max(copter_x)+1.1*pL]
    lim_y = [min(copter_y)-1.1*pL max(copter_y)+1.1*pL]
    lim_z = [min(copter_z)-1.1*pL max(copter_z)+1.1*pL]
    axis([ lim_x lim_y lim_z]);

    axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on;
end
% 
% lim_x = [min(eta(1,:))-0.5 max(eta(1,:))+0.5]
% lim_y = [min(eta(2,:))-0.5 max(eta(2,:))+0.5]
% lim_z = [min( [min(eta(3,:))-0.1 min(eta_l(3,:))]) max([max(eta(3,:)+0.1) max(eta_l(3,:))])] 
%
% 
% axis([ lim_x lim_y lim_z]);
% Plot the result of convolution

figure(2); clf
hold on;
plot(time, log.get('r_preconv', 1), 'b');
plot(time, log.get('r', 1), 'b--');


figure(1);
mainFHandle = figure(1);
% Plot nicely

old_copter = [];
old_plot = [];
old_load = [];
skips = 10;
timeShifter = 1;

t = 1;

if doOnlyLastFrame
    t = length(time)-1
end

movie = 0;

if movie
   M =  cell(1, length(time));
end

colors = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

colors = mat2cell(colors, ones(1, size(colors, 1)));
while t < length(time)
    
    %profile on
    timeToDraw = tic;
  
    p.deleteOld();
    
    % PLOT MAIN COPTER VIEW
    [on, subcell] = p.isActive('main');
    if on
        subplot(subcell{1}{:})
        

        if all(ishandle(old_copter))
            delete(old_copter);
        end
        
        % Plot current copter position

        
        p.plot3(log.get('Eta', 1, t), log.get('Eta', 2, t), log.get('Eta', 3, t), 'ro');
        
        % Plot line to the load
        x_copter = log.get('Eta', 1, 1:t);
        y_copter = log.get('Eta', 2, 1:t);
        z_copter = log.get('Eta', 3, 1:t);
        
        
        x_load = x_copter + pL*sin(-log.get('Eta', 3, 1:t));
        y_load = y_copter + pL*cos(-log.get('Eta', 3, 1:t));
        y_load = y_copter + pL*cos(-log.get('Eta', 3, 1:t));
        
        
        phi_L_now = log.get('Eta', 4, t);
        theta_L_now = log.get('Eta', 5, t);
        
        p_load = log.get('Eta', 1:3, t) + f_5dof_Rload(phi_L_now, theta_L_now)*[0; 0; pL];

        h_line = p.line([x_copter(end) p_load(1)], [y_copter(end) p_load(2)], [z_copter(end) p_load(3)]);
        set(h_line, 'color', 'k');
        
        % Plot load
        p.plot3(p_load(1), p_load(2), p_load(3), 'go');
        
        % Plot history of load
        %p.plot(x_load, y_load, 'k--');
        
        % Plot history copter
        p.plot3(x_copter, y_copter, z_copter, 'b--');
        
        % Plot desired position
        p.plot3(r(1,1:t), r(2,1:t), r(3,1:t), 'g');

        

        axis equal;
        axis([ lim_x lim_y lim_z]);
        if movie
            % Trhee moving
            %axis([-2.0913    1.0328   -1.5262    1.2083   -1.8269    1.3284]);
            
            % two demo
            %axis([-3.1999    1.3166   -0.8530    0.9999   -1.7931    1.6479]);
            
            
            % Wierd one no control
            %axis([ -2.8706    1.2828   -4.0354    2.8947   -1.7311    6.2465]);
            
            % 3 tracking test
            %axis([    1.0000   12.8226    0.6557    3.3443    1.0010    3.8016]);

            
            
            
        end
    end
    
    
    % PLOT Theta_L and phi_L
    [on, subcell] = p.isActive('theta');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;
        
        hold on;
        axis fill
        xlim([0 tend]);
        title('Theta_L and Phi_L'); xlabel('time [s]'); ylabel('theta_L [deg]');
        a = p.plot(time(1:t), pi\180*log.get('Eta', 4, 1:t));
        %a.Color = colors{1};
        p.plot(time(1:t), pi\180*log.get('Eta', 5, 1:t));
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
    end
    
    % PLOT pos_error
    [on, subcell] = p.isActive('pos_error');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;
        
        hold on;
        axis fill
        xlim([0 tend])
        title('Position error'); xlabel('time [s]'); ylabel('error');
        p.plot(time(1:t), log.get('Eta', 1,1:t) - r(1,1:t));
        p.plot(time(1:t), log.get('Eta', 2,1:t) - r(2,1:t));
        p.plot(time(1:t), log.get('Eta', 3,1:t) - r(3,1:t));
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
        grid on;
        
    end
    
    % PLOT pos_error
    [on, subcell] = p.isActive('vel_error');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;
        
        hold on;
        axis fill
        xlim([0 tend])
        title('Velocity error'); xlabel('time [s]'); ylabel('error');
        p.plot(time(1:t), log.get('Nu', 1,1:t) - dr(1,1:t));
        p.plot(time(1:t), log.get('Nu', 2,1:t) - dr(2,1:t));
        p.plot(time(1:t), log.get('Nu', 3,1:t) - dr(3,1:t));
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
        grid on;
        
    end
    
    % PLOT vel norm
    [on, subcell] = p.isActive('vel');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;
        
        hold on;
        axis fill
        xlim([0 tend])
        title('Velocity norm'); xlabel('time [s]'); ylabel('norm of velocity');
        p.plot(time(1:t), sqrt(sum(abs(log.get('Nu', 1:3,1:t)).^2,1)));
        p.plot(time(1:t), sqrt(sum(abs(log.get('d2r_postconv', ':',1:t)).^2,1)));
        
        p.plot(time(1:t), log.get('feedback_active', 1, 1:t));
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
        grid on;
        if t == 1
            legend('Vel Norm', 'accref', 'delayedGain'); 
        end
        
    end
    
    % PLOT pos_error
    [on, subcell] = p.isActive('acc_error');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;
        
        hold on;
        axis fill
        xlim([0 tend])
        title('Acc error'); xlabel('time [s]'); ylabel('error');
        p.plot(time(1:t), log.get('dNu', 1,1:t) - d2r(1,1:t));
        p.plot(time(1:t), log.get('dNu', 2,1:t) - d2r(2,1:t));
        p.plot(time(1:t), log.get('dNu', 3,1:t) - d2r(3,1:t));
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
        grid on;
        
    end
    
    [on, subcell] = p.isActive('alpha_danglular');
    if on
       subplot(subcell{1}{:})
       
       a = gca; a.ColorOrderIndex = 1;
       
       hold on;
       axis fill
       xlim([0 time(end)]);
       title('Alpha vs angular rates. ');
       
       
       p.plot(time(1:t), log.get('Nu', 4, 1:t), '-');
       a.ColorOrderIndex = 2;
       p.plot(time(1:t), log.get('alpha_45', 1, 1:t), '--');
       
       
       p.plot(time(1:t), log.get('Nu', 5, 1:t), '-');
       a.ColorOrderIndex = 4;
       p.plot(time(1:t), log.get('alpha_45', 2, 1:t), '--');
       
    end
    
        % PLOT pos_error
    [on, subcell] = p.isActive('loadEnergy');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;
        
        hold on;
        axis fill
        xlim([0 tend])
        title('Load Energy'); xlabel('time [s]'); ylabel('energy');
        %p.plot(time(1:t), loadEnergy(1,1:t));
        p.plot(time(1:t), loadEnergy(2,1:t));
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
        grid on;
        
    end
    
    % PLOT taur
    [on, subcell] = p.isActive('tau');
    if on
        subplot(subcell{1}{:})
        
        hold on;
        axis fill
        xlim([0 tend])
        title('Control input'); xlabel('time [s]'); ylabel('tau');
        p.plot(time(1:t), log.get('tau', 1,1:t), 'b');
        p.plot(time(1:t), log.get('tau', 2,1:t), 'g');
        p.plot(time(1:t), log.get('tau', 3,1:t), 'r');
        p.plot(time(1:t), log.get('tau_45', 1,1:t), 'k--');
        p.plot(time(1:t), log.get('tau_45', 2,1:t), 'k--');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
    end
    
    % PLOT m_l_hat
    [on, subcell] = p.isActive('m_L_hat');
    if on
        subplot(subcell{1}{:})
        
        hold on;
        axis fill
        xlim([0 tend]);
        title('Payload mass estimate'); xlabel('time [s]'); ylabel('m_L_hat');
        p.plot(time(1:t), log.get('m_L_hat', 1, 1:t));
        p.plot(time, pm_L*ones(size(time)), '--');
        
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
    end
    
    
    % PLOT integral action
    [on, subcell] = p.isActive('z2_integral');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;

        
        hold on;
        axis fill
        xlim([0 tend]);
        title('Integral action'); xlabel('time [s]'); ylabel('integral');
        p.plot(time(1:t), log.get('z2_integral', 1, 1:t));
        p.plot(time(1:t), log.get('z2_integral', 2, 1:t));
        p.plot(time(1:t), log.get('z2_integral', 3, 1:t));
        
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
    end
    
    
    
    % PLOT position
    [on, subcell] = p.isActive('pos');
    if on
        subplot(subcell{1}{:})
        a = gca; a.ColorOrderIndex = 1;
        
        hold on;
        axis fill
        xlim([0 tend])
        title('Position'); xlabel('time [s]'); ylabel('error');
        p.plot(time(1:t), log.get('Eta', 1,1:t), 'b');
        p.plot(time(1:t), r(1,1:t), 'b--');
        p.plot(time(1:t), log.get('Eta', 2,1:t), 'r');
        p.plot(time(1:t), r(2,1:t), 'r--');
        p.plot(time(1:t), log.get('Eta', 3,1:t), 'g');
        p.plot(time(1:t), r(3,1:t), 'g--');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_r_storage(4, 1:t), 'r');
        %p.plot(time(1:t), pi\180*sys.Copters(infoCopter).Eta_comp_storage(4, 1:t), 'g');
        
        grid on;
        
    end
    
    infoCopter = 1;
    
  
    
    drawnow;
    %profile viewer
    %return
    
    if movie
        M{t} = getframe(gcf);
    end
    
    timeToDraw = toc(timeToDraw);
    
    steps = int32(timeToDraw/h);
    if steps < 1
        steps = 1;
    end
    %steps = 1;
    
    if movie
        steps = 5;
    end
    
    t = t + steps;
    
    %pause(timeShifter*skips*h)
    if timeToDraw < h
        pause(h-timeToDraw)
    end
    
end

end % interactive plot
%%
% Plot some of the controller effects
% figure(3); clf
% hold on;
% plot(time, log.get('tau_beta_1'), '--')
% plot(time, log.get('tau_beta_2'), '.-')
% plot(time, log.get('tau_main'))
