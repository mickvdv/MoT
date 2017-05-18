format long g

% simulation parameters
N = 10;

% load basket from excel file

% environmental variables
q = 0.0275;
sig = 0.3; % per year
T = 10;

% product specifications
premium = 150;
global payments_per_year;
payments_per_year = 12;
participation_rate = 0.75;
cap_rate = 1.15; % on yearly basis

% initialization
S0 = 100;
StockPath = zeros(N, T*payments_per_year+1);
StockPath(:, 1) = S0;

% simulate stock path
for sim_i = 1:N
    for month_i = 2:(T*payments_per_year)+1
%         This is the simulation of the next stock price
        r = RiskFreeRateInterpolation(10);
        StockPath(sim_i, month_i) = normrnd(StockPath(sim_i, month_i-1)*exp(r*(1/payments_per_year)), sig);
    end
end

% simulate the call prices
CallPrices = zeros(N, T*payments_per_year);
ShortStrikes = zeros(N, T*payments_per_year);
LongStrikes = zeros(N, T*payments_per_year);
CallAmounts = zeros(N, T*payments_per_year);

for sim_i = 1:N
    for month_i = 1:(T*payments_per_year)
        t = T - ((month_i-1) * (1/payments_per_year));
        r = RiskFreeRateInterpolation(t);
        
        amount_of_stock = (premium * participation_rate) / StockPath(sim_i, month_i);
        CallAmounts(sim_i, month_i) = amount_of_stock;
        
        K_long = StockPath(sim_i, month_i);
        call_long = bsm_call(r, q, StockPath(sim_i, month_i), StockPath(sim_i, month_i), t, sig);
        LongStrikes(sim_i, month_i) = K_long;
        
        K_short = StockPath(sim_i, month_i) +  ((cap_rate - 1) * t) * StockPath(sim_i, month_i);
        ShortStrikes(sim_i, month_i) = K_short;
        call_short = bsm_call(r, q, StockPath(sim_i, month_i), K_short, t, sig);
        
        
        CallPrices(sim_i, month_i) = call_long * amount_of_stock - call_short * amount_of_stock;
    end
end

% simulate the bond prices
BondPrices = zeros(N, T*payments_per_year);
BondFaceValues = zeros(N, T*payments_per_year);

for sim_i = 1:N
    for month_i = 1:(T*payments_per_year)
        t = T - ((month_i-1) * (1/payments_per_year));
        r = RiskFreeRateInterpolation(t);
        
%         use whatever money is left to buy the bond
        left_over = premium - CallPrices(sim_i, month_i);
        BondFaceValues(sim_i, month_i) = left_over * exp(r * t);
        BondPrices(sim_i, month_i) = BondFaceValues(sim_i, month_i) * exp(-r * t);
        
    end
end


% Now we have the products bought, we can calculate the payoff
PayOff = zeros(N, 1);

for sim_i = 1:N
    PayOff(sim_i) = sum(BondFaceValues(sim_i, :),2);
    final_value_of_stock = StockPath(sim_i, T*payments_per_year+1);
    for month_i = 1:(T*payments_per_year)
        PayOff(sim_i) = PayOff(sim_i) + CallAmounts(sim_i, month_i) * (max(0, final_value_of_stock - LongStrikes(sim_i, month_i)) - max(0, final_value_of_stock - ShortStrikes(sim_i, month_i)));
    end
end

% Reporting
PortfolioPrices = BondPrices + CallPrices;
% disp(['Total price of the porfolio ' num2str(sum(PortfolioPrices, 2))]);
% disp(['Principal that is protected ' num2str(sum(BondFaceValues, 2))]);
% disp(['Customer premium ' num2str(payments_per_year * T * premium)]);
disp(['Protection rate: ' num2str(sum(BondFaceValues(1,:),2)/(payments_per_year * T * premium))]);


% plot(risk_free_rate);

function StockPath(S0, T, mu, sig)
   
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
    global payments_per_year;
    % Risk free rate interpolation
    risk_free_rate = zeros(1, 20 * payments_per_year);
    risk_free_rate(:,:) = NaN;
    risk_free_rate(20*payments_per_year) = 2.71;
    risk_free_rate(10*payments_per_year) = 2.33;
    risk_free_rate(7*payments_per_year) = 2.13;
    risk_free_rate(5*payments_per_year) = 1.84;
    risk_free_rate(3*payments_per_year) = 1.48;
    risk_free_rate(2*payments_per_year) = 1.28;
    risk_free_rate(1*payments_per_year) = 1.09;
    risk_free_rate(0.5*payments_per_year) = 0.98;
    risk_free_rate(0.25*payments_per_year) = 0.83;
    risk_free_rate((1/12)*payments_per_year) = 0.67;

    idxValid = ~isnan(risk_free_rate);

    x = linspace(1,20*12,20*12);
    p = polyfit(x(idxValid),risk_free_rate(idxValid),3);

    risk_free_rate_interpolation = polyval(p,x);
    t_in_months = round( t * payments_per_year);
    r = risk_free_rate_interpolation(t_in_months)/100;
end



