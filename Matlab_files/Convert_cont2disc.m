function [A_d, B_d] = Convert_cont2disc(A,B,h)
%[A_d, B_d] = Convert_cont2disc(A,B,h)
%Convert a continuous time system to a discrete: 
%x_dot = A_c*x + B_c*u to x_(k+1) = A_d*x(k) + B_d*u(k)
%
%Author: Pål Mathisen - 18.10.2017

    A_d = expm(A*h);
    conv_phi = eye(length(A));
    A_i = A;
    for i = 1:100
        conv_phi = conv_phi + (A_i)/factorial(i+1)*(h^i);
        A_i = A_i * A;
    end
    B_d = h*conv_phi*B;
end