function [y] = h_EK(x, sig)
%[y] = h_EK(x, sig)
    N_measurment            = size(sig,2); 
    y                       = zeros(N_measurment,1);
    for i = 1:N_measurment
        y(i)                = norm(x(1:3)-sig(:,i));
    end
    
end