function [P_pri_next, x_hat_next, x_hat_post, error] = EKF_range(x_hat, sig, y, u, P_pri, A, B, R, Q, H_EK)
    %[P_pri_next, x_hat_next, x_hat_post] = EKF_range(x_hat, sig, y, u, P_pri, A, B, R, Q, H_EK)
    
    %Prework
    n = size(x_hat,1);
    y_hat = h_EK(x_hat, sig);
    H                       = H_EK(x_hat(1:3), sig, y_hat);
    if rank(obsv(A,H)) ~= 9
        error('Not observable');
    end
    
    %Corrector
    K                       = P_pri*H'/(H*P_pri*H' + R);
    error                   = y-y_hat;
    x_hat_post              = x_hat + K*(error);
    P_post                  = (eye(n) - K*H)*P_pri*(eye(n) - K*H)' + K*R*K';
    P_post                  = (P_post + P_post')/2;
    
    % Predictor
    x_hat_next              = A*x_hat_post + B*u;
    P_pri_next              = A*P_post*A' + Q;
end