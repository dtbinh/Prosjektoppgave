function [P_pri_next, x_hat_next, x_hat_post, K] = KF_range(x_hat, x_m, P_pri, A, R, Q, H)
    %[P_pri_next, x_hat_next, x_hat_post, K] = KF_range(x_hat, x_m, P_pri, A, R, Q, H)

    %Corrector
    K                       = P_pri*H'/(H*P_pri*H' + R);
    x_hat_post              = x_hat + K*(x_m-x_hat);
    P_post                  = (eye(6) - K*H)*P_pri*(eye(6) - K*H)' + K*R*K';
    
    % Predictor
    x_hat_next              = A*x_hat_post;
    P_pri_next              = A*P_post*A' + Q;
end