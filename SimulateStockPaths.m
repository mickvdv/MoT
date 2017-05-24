function StockPaths = SimulateStockPaths(S0, T, dT, mu, sig, N)
%    http://www.investopedia.com/articles/07/montecarlo.asp
    t_size = round(T/dT);
    StockPaths = zeros(N, t_size+1);
    StockPaths(:, 1) = S0;
    for sim_i = 1:N
        StockPaths(sim_i, :) = GenerateStockPath(S0, T, dT, mu, sig);
    end           
end