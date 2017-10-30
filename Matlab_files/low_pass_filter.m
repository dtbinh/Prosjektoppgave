function lp_sig = low_pass_filter(sig,T)
    lp_sig      = zeros(size(sig, 1), size(sig, 2));
    lp_sig(:,1) = sig(:,1);
    len         = size(sig,2);
    for k = 2:len
        lp_sig(:,k) = lp_sig(:,k-1) + 1/T*(sig(:,k) - lp_sig(:,k-1) );
    end
end
