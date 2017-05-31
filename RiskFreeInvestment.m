

final_value = 0;
load('rfr.mat', 'risk_free_rate_interpolation');

for T_to_mat = 1:(1/12):10
    r = RiskFreeRateInterpolation(T_to_mat);
    final_value = final_value + 150 * exp(r*T_to_mat);
end

final_value / (150*120)