function [delta, gamma, theta, vega, rho] = bsm_call_greeks(r, q, S0, K, T, sig)
    S0 = S0*exp(-q * T);
    d1 = (log(S0/K)+(r+(sig^2/2))*T)/(sig*sqrt(T));
    d2 = d1 - sig*sqrt(T);
    
    n_pdf_d1 = normpdf(d1);
    n_pdf_d2 = normpdf(d2);
    n_cdf_d2 = normcdf(d2);
    
    theta = - (S0*n_pdf_d1*sig)/(2*sqrt(T))-r*K*exp(-r*T)*n_cdf_d2;
    delta = n_pdf_d1;
    gamma = n_pdf_d1/(S0*sig*sqrt(T));
    vega = S0*sqrt(T)*n_pdf_d1;
    rho = K*T*exp(-r*T)*n_pdf_d2;
end
