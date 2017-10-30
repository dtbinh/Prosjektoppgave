function [ states_, P_, y_tilde ] = step_gnss_uwb( Ts, s_states, P_, c_gnss_meas, c_uwb_meas, Q, nav_frame )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % constant
    I3 = eye(3);
    I2 = eye(2);
    I8 = eye(8);
    z3 = zeros(3);
    z32 = zeros(3,2);
    z23 = zeros(2,3);
    Ad = [  I3  Ts*I3   z32; 
            z3  I3      z32;
            z23 z23     I2];
    Bd = [  Ts^2/2*I3   z32; 
            Ts*I3       z32;
            z23         Ts*I2];
        
    % input states
    x_ = [...
        s_states.pos; 
        s_states.velo;
        s_states.beta_gnss;
        s_states.beta_uwb;
        ];
    R_en = s_states.R_en;
    
    %% meas corrections
    y = [];
    y_hat = [];
    R = [];
    C0 = [];
    C0_NED = [];
    p_eb_e_pred = x_(1:3);
    
    % GNSS 
    beta_gnss = s_states.beta_gnss;
    N = length(c_gnss_meas);
    for k = 1:N
        s_gnss = c_gnss_meas{k};
        p_eb_e_k = s_gnss.p_eb_e_k;
        y_k = s_gnss.y;
        r_k = s_gnss.std^2;
        rho_hat_k = norm( p_eb_e_pred - p_eb_e_k);
        LOS_k = (p_eb_e_pred - p_eb_e_k )/rho_hat_k;
        if strcmp(nav_frame,'ned')
            LOS_NED = LOS_k;
            C0_NED = [C0_NED; LOS_k' 1 0];
            C0 = [C0; LOS_k'];
            y = [y; y_k];
            y_hat = [y_hat; rho_hat_k + beta_gnss];
            R = blkdiag(R, r_k);
        elseif strcmp(nav_frame,'ecef')
            LOS_NED = R_en'*LOS_k;
            elevation_angle = asind( LOS_NED(3) );
            if elevation_angle >= 15
                C0_NED = [C0_NED; LOS_NED' 1 0];
                C0 = [C0; LOS_k'];
                y = [y; y_k];
                y_hat = [y_hat; rho_hat_k + beta_gnss];
                R = blkdiag(R, r_k);
            end
        else
            error('Wrong nav frame %s', nav_frame)
        end
    end
    m_gnss = size(C0, 1);
    
    % uwb 
    beta_uwb = s_states.beta_uwb;
    N = length(c_uwb_meas);
    for k = 1:N
        s_uwb = c_uwb_meas{k};
        p_eb_e_k = s_uwb.p_eb_e_k;
        y_k = s_uwb.y;
        r_k = s_uwb.std^2;
        rho_hat_k = norm( p_eb_e_pred - p_eb_e_k);
        LOS_k = (p_eb_e_pred - p_eb_e_k )/rho_hat_k;
        if strcmp(nav_frame,'ned')
            LOS_NED = LOS_k;
            C0_NED = [C0_NED; LOS_k' 0 1];
        elseif strcmp(nav_frame,'ecef')
            LOS_NED = R_en'*LOS_k;
            C0_NED = [C0_NED; LOS_NED' 0 1];
        else
            error('Wrong nav frame %s', nav_frame)
        end
        C0 = [C0; LOS_k'];
        y = [y; y_k];
        y_hat = [y_hat; (rho_hat_k + beta_uwb)];
        R = blkdiag(R, r_k);
    end
    m_uwb = size(C0, 1)-m_gnss;
    
    C = [C0 zeros( m_gnss+m_uwb, 3) [ones(m_gnss, 1); zeros(m_uwb, 1)] [zeros(m_gnss, 1); ones(m_uwb, 1)]  ];
    K = P_*C'/(C*P_*C' + R);
    y_tilde = y - y_hat;
    x_hat = x_ + K*y_tilde;
    P_hat = (I8- K*C)*P_*(I8- K*C)' + K*R*K';
    P_hat = (P_hat + P_hat')/2;
    % predictions
    x_ = Ad*x_hat;
    Qd = Q*sqrt(1/Ts);
    P_ = Ad*P_hat*Ad' + Bd*Qd*Bd';
    P_ = (P_ + P_')/2;
    
    % output    
    states_.pos = x_(1:3); 
    states_.velo = x_(4:6);
    states_.beta_gnss = x_(7);
    states_.beta_uwb = x_(8);
    if strcmp(nav_frame,'ned')
        lat_hat = 0; lon_hat = 0; height_hat = 0;
    elseif strcmp(nav_frame,'ecef')
        [lat_hat, lon_hat, height_hat] = ecef2lla( x_(1), x_(2), x_(3) );
    end
    states_.lat = lat_hat;
    states_.lon = lon_hat;
    states_.height = height_hat; 
    if strcmp(nav_frame,'ned')
        states_.R_en = I3;
    elseif strcmp(nav_frame,'ecef')
        states_.R_en = Rll(lon_hat,lat_hat);
    end
    if isempty( c_uwb_meas ) && isempty( c_gnss_meas )
        error('not measurement. DOP calc to possible')
    elseif isempty( c_uwb_meas ) 
        C0_NED = C0_NED(:,1:end-1);
    elseif isempty( c_gnss_meas )
        C0_NED = [ C0_NED(:,1:end-2) C0_NED(:,end) ];
    end   
    G_DOP = inv(C0_NED'*C0_NED);
    states_.HDOP = sqrt( G_DOP(1,1) + G_DOP(2,2) );
    states_.VDOP = sqrt( G_DOP(3,3) );
    if ~isempty( c_gnss_meas ) 
        states_.TDOP_GNSS = sqrt( G_DOP(4,4) );
    else
        states_.TDOP_GNSS = 0;
    end
    if ~isempty( c_gnss_meas )  && ~isempty( c_uwb_meas ) 
        states_.TDOP_UWB = sqrt( G_DOP(5,5) );
    elseif isempty( c_gnss_meas )  && ~isempty( c_uwb_meas ) 
        states_.TDOP_UWB = sqrt( G_DOP(4,4) );
    else
        states_.TDOP_UWB = 0;
    end
end

