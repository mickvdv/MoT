function c = bsm_call(r, q, S0, K, T, sig)
    S0 = S0*exp(-q * T);
    d1 = (log(S0/K)+(r+(sig^2/2))*T)/(sig*sqrt(T));
    d2 = d1 - sig*sqrt(T);
    c = S0 * normcdf(d1) - K * exp(-r*T)*normcdf(d2);
%     p = K * exp(-r*T)*normcdf(-d2)-S0*normcdf(-d1)
%     dc = normcdf(d1);
%     dp = normcdf(d1)-1;
%     g = normpdf(d1)/(S0*sig*sqrt(T));
end
