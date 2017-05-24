function r = RiskFreeRateInterpolation(t)
     global risk_free_rate_interpolation;
     global payments_per_year;
    t_in_months = round( t * payments_per_year);
    r = risk_free_rate_interpolation(t_in_months);
end
