format long g

% simulation parameters
N = 50;

% for loading the basket
filename = 'basket.xlsx';

% Load the most recent risk free interpolation
load('rfr.mat', 'risk_free_rate_interpolation');

% product specifications
premium = 150;
global payments_per_year;
payments_per_year = 12;
T = 10;
dT = 1/payments_per_year;

% load the basket
% Array contents: (S0, sig, null, null, q)
basket = xlsread(filename)
basket_size = size(basket, 1);

CapRates = 1:0.025:1.5;
ParticipationRates = 0.2:0.025:1.2;


ProtectionRates = zeros(size(CapRates,2), size(ParticipationRates, 2));

for cap_i = 1:size(CapRates,2)
    cap_rate = CapRates(1, cap_i);
    for participation_i = 1:size(ParticipationRates, 2)
        participation_rate = ParticipationRates(participation_i);
        % calculate the call prices
        CallPrices = zeros(basket_size, T*payments_per_year);
        CallAmounts = zeros(basket_size, T*payments_per_year);

        for stock_i = 1 : basket_size
        %     for each stock
            S0 = basket(stock_i, 1);
            sig = basket(stock_i, 2);
            q = basket(stock_i, 5);

            for month_i = 1:(T/dT)
                sim_T = T - (month_i - 1) * dT;
                r = RiskFreeRateInterpolation(sim_T);
                StockPaths = SimulateStockPaths(S0, sim_T, dT, r - q, sig, N);

        %         get the expected values of the calls
                ExpectedLongCallValue = ExpectedCallValueFromStockPaths(StockPaths, S0, sim_T, dT, N);
                K_short = S0 +  ((cap_rate - 1) * sim_T) * S0;
                ExpectedShortCallValue = ExpectedCallValueFromStockPaths(StockPaths, K_short, sim_T, dT, N);

                amount_of_stock = (1/basket_size) * (premium * participation_rate) / S0;

        %         save prices and amounts
                CallPrices(stock_i, month_i) = ExpectedLongCallValue * amount_of_stock - ExpectedShortCallValue * amount_of_stock;
                CallAmounts(stock_i, month_i) = amount_of_stock;

        %         prepare next simulation step
                S0 = ExpectedValueFromStockPaths(StockPaths, dT, dT, N);     
            end
        end

        % Calculate the bonds
        BondPrices = 150 - sum(CallPrices, 1);
        BondFaceValues = zeros(1, (T/dT));

        for month_i = 1:(T/dT)
            sim_T = T - (month_i - 1) * dT;
            r = RiskFreeRateInterpolation(sim_T);
            BondFaceValues(1, month_i) = BondPrices(1, month_i) * exp(r * sim_T);
        end

        % Protection
        Protection = sum(BondFaceValues, 2);
        ProtectionRates(cap_i, participation_i) = Protection / (premium * payments_per_year * T);
    end
end


function e = ExpectedCallValueFromStockPaths(StockPaths, K, forwardT, dT, N)
    t_index = round(forwardT/dT);
    r = RiskFreeRateInterpolation(forwardT);
    call_values = exp(-r * forwardT) * max(StockPaths(:, t_index) - K, 0);
    e = sum(call_values)/N;
end

function e = ExpectedValueFromStockPaths(StockPaths, forwardT, dT, N)
    t_index = round(forwardT/dT) + 1;
    col_sum = sum(StockPaths, 1);
    e = col_sum(t_index)/N;
end
    

function StockPaths = SimulateStockPaths(S0, T, dT, mu, sig, N)
%    http://www.investopedia.com/articles/07/montecarlo.asp
    t_size = round(T/dT);
    StockPaths = zeros(N, t_size);
    StockPaths(:, 1) = S0;
    for sim_i = 1:N
        for step_i = 2:(T/dT) + 1 %(plus 1 is for the last month)
            StockPaths(sim_i, step_i) = max(StockPaths(sim_i, step_i - 1) + (StockPaths(sim_i, step_i - 1) * (mu * dT + sig * randn() * sqrt(dT))), 0);
        end
    end           
end


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


function r = RiskFreeRateInterpolation(t)
    global risk_free_rate_interpolation;
    global payments_per_year;
    t_in_months = round( t * payments_per_year);
    r = risk_free_rate_interpolation(t_in_months);
end



