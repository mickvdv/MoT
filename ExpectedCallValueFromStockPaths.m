function e = ExpectedCallValueFromStockPaths(StockPaths, K, forwardT, dT, N)
    t_index = round(forwardT/dT);
    r = RiskFreeRateInterpolation(forwardT);
    call_values = exp(-r * forwardT) * max(StockPaths(:, t_index) - K, 0);
    e = sum(call_values)/N;
end