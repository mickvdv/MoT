function e = ExpectedValueFromStockPaths(StockPaths, forwardT, dT, N)
    t_index = round(forwardT/dT) + 1;
    col_sum = sum(StockPaths, 1);
    e = col_sum(t_index)/N;
end