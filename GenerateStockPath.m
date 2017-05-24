function StockPath = GenerateStockPath(S0, T, dT, mu, sig)
    StockPath = zeros((T/dT));
    StockPath(1) = S0;
    for step_i = 2:(T/dT) + 1 %(plus 1 is for the last month)
        StockPath(step_i) = max(StockPath(step_i - 1) + (StockPath(step_i - 1) * (mu * dT + sig * randn() * sqrt(dT))), 0);
    end
end