format long g

% simulation parameters
N = 100;
% for loading the basket
filename = 'basket.xlsx';

% Load the most recent risk free interpolation
load('rfr.mat', 'risk_free_rate_interpolation');
global risk_free_rate_interpolation;

% product specifications
premium = 150;
global payments_per_year;
payments_per_year = 12;
participation_rate = 0.875;
cap_rate = 1.075; % on yearly basis
T = 10;
dT = 1/payments_per_year;

% load the basket
% Array contents: (S0, sig, null, null, q)
basket = xlsread(filename)
basket_size = size(basket, 1);

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
ProtectionRate = Protection / (premium * payments_per_year * T)

% CallPrices



% 
% 
% % simulate the call prices
% CallPrices = zeros(N, T*payments_per_year);
% ShortStrikes = zeros(N, T*payments_per_year);
% LongStrikes = zeros(N, T*payments_per_year);
% CallAmounts = zeros(N, T*payments_per_year);
% 
% for sim_i = 1:N
%     for month_i = 1:(T*payments_per_year)
%         t = T - ((month_i-1) * (1/payments_per_year));
%         r = RiskFreeRateInterpolation(t);
%         
%         amount_of_stock = (premium * participation_rate) / StockPath(sim_i, month_i);
%         CallAmounts(sim_i, month_i) = amount_of_stock;
%         
%         K_long = StockPath(sim_i, month_i);
%         call_long = bsm_call(r, q, StockPath(sim_i, month_i), StockPath(sim_i, month_i), t, sig);
%         LongStrikes(sim_i, month_i) = K_long;
%         
%         K_short = StockPath(sim_i, month_i) +  ((cap_rate - 1) * t) * StockPath(sim_i, month_i);
%         ShortStrikes(sim_i, month_i) = K_short;
%         call_short = bsm_call(r, q, StockPath(sim_i, month_i), K_short, t, sig);
%         
%         
%         CallPrices(sim_i, month_i) = call_long * amount_of_stock - call_short * amount_of_stock;
%     end
% end
% 
% % simulate the bond prices
% BondPrices = zeros(N, T*payments_per_year);
% BondFaceValues = zeros(N, T*payments_per_year);
% 
% for sim_i = 1:N
%     for month_i = 1:(T*payments_per_year)
%         t = T - ((month_i-1) * (1/payments_per_year));
%         r = RiskFreeRateInterpolation(t);
%         
% %         use whatever money is left to buy the bond
%         left_over = premium - CallPrices(sim_i, month_i);
%         BondFaceValues(sim_i, month_i) = left_over * exp(r * t);
%         BondPrices(sim_i, month_i) = BondFaceValues(sim_i, month_i) * exp(-r * t);
%         
%     end
% end
% 
% 
% % Now we have the products bought, we can calculate the payoff
% PayOff = zeros(N, 1);
% 
% for sim_i = 1:N
%     PayOff(sim_i) = sum(BondFaceValues(sim_i, :),2);
%     final_value_of_stock = StockPath(sim_i, T*payments_per_year+1);
%     for month_i = 1:(T*payments_per_year)
%         PayOff(sim_i) = PayOff(sim_i) + CallAmounts(sim_i, month_i) * (max(0, final_value_of_stock - LongStrikes(sim_i, month_i)) - max(0, final_value_of_stock - ShortStrikes(sim_i, month_i)));
%     end
% end
% 
% % Reporting
% PortfolioPrices = BondPrices + CallPrices;
% % disp(['Total price of the porfolio ' num2str(sum(PortfolioPrices, 2))]);
% % disp(['Principal that is protected ' num2str(sum(BondFaceValues, 2))]);
% % disp(['Customer premium ' num2str(payments_per_year * T * premium)]);
% disp(['Protection rate: ' num2str(sum(BondFaceValues(1,:),2)/(payments_per_year * T * premium))]);
% 
% 
% % plot(risk_free_rate);

function e = ExpectedCallValueFromStockPaths(StockPaths, K, forwardT, dT, N)
    t_index = round(forwardT/dT);
    call_values = max(StockPaths(:, t_index) - K, 0);
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
            StockPaths(sim_i, step_i) = max(StockPaths(sim_i, step_i - 1) + (StockPaths(sim_i, step_i - 1) * (mu * dT + sig * normrnd(0,1) * sqrt(dT))), 0);
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



