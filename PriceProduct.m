% product specifications
premium = 150;
global payments_per_year;
payments_per_year = 12;

global T;
T = 10;
dT = 1/payments_per_year;

cap_rate = 1.10;
participation_rate = 0.95;

% Load bond face values
load('bonds.mat', 'BondFaceValues');
load('rfr.mat', 'risk_free_rate_interpolation');

% basket parameters
global sig;
global q;
S0 = 68.63;
sig = 0.1567;
q = 0.0203;

% pricing
global r_stock;
r_stock = 0.00822600659112768;
month_duration = 21; % trading days per month


StockPath = GenerateStockPath(S0, T, (1/(month_duration*payments_per_year)), (r_stock-q), sig);

LongCallStrikes = zeros(1, T/dT);
LongCallAmounts = zeros(1, T/dT);
ShortCallAmounts = zeros(1, T/dT);
ShortCallStrikes = zeros(1, T/dT);
Deltas = zeros(1, T/dT);
Gammas = zeros(1, T/dT);
Thetas = zeros(1, T/dT);
Vegas = zeros(1, T/dT);
Rhos = zeros(1, T/dT);

for month_i = 1:(T/dT)
    day_i = ((month_i-1) * month_duration) + 1;
    sim_T = T - (month_i - 1) * dT;
    S0 = StockPath(1, day_i);
    amount_of_stock = (premium * participation_rate) / S0;
    
    LongCallAmounts(1, month_i) = amount_of_stock;
    ShortCallAmounts(1, month_i) = amount_of_stock;
    LongCallStrikes(1, month_i) = S0;
    ShortCallStrikes(1, month_i) = S0 * cap_rate ^ (sim_T);    
end

values = zeros(1, 10 * month_duration*payments_per_year);
for day_i = 1:(size(StockPath, 2)-1)
%     day_i
    S0 = StockPath(1, day_i);
    
     [present_value, Delta, Gamma, Theta, Vega, Rho] = CalculateProductPrice(S0, LongCallStrikes, LongCallAmounts, ShortCallStrikes, ShortCallAmounts, BondFaceValues, T - (day_i/(month_duration*payments_per_year)));
     values(1, day_i) = present_value;
     Deltas(1, day_i) = Delta;
     Gammas(1, day_i) = Gamma;
     Thetas(1, day_i) = Theta;
     Vegas(1, day_i) = Vega;
     Rhos(1, day_i) = Rho;
end

figure
subplot(3,1,1)       % add first plot in 2 x 1 grid
plot(StockPath);
yyaxis right
xlabel('Time (in days)')
plot(values);

subplot(3,1,2)       % add second plot in 2 x 1 grid
hold on
plot(Deltas)
plot(Gammas)
hold off
xlabel('Time (in days)')
legend('Delta', 'Gamma');

subplot(3,1,3)
hold on
plot(Thetas)
plot(Vegas)
plot(Rhos)
hold off
xlabel('Time (in days)')
legend('Theta', 'Vega', 'Rho');


PayOff = values(1, 252*10) / sum(BondFaceValues)

function [present_value, Delta, Gamma, Theta, Vega, Rho] = CalculateProductPrice(S0, LongCallStrikes, LongCallAmounts, ShortCallStrikes, ShortCallAmounts, BondFaceValues, sim_T)
%     T in year, till maturity
    global payments_per_year;
    global r_stock;
    global q;
    global sig;
    global T;
    
    
    dT = 1/payments_per_year;
%     In which month are we currently
    current_month = min(floor((T-sim_T) * payments_per_year)+1, T * payments_per_year);
    OptionValue = 0;
    BondValue = 0;
%     bsm_call(r, q, S0, K, T, sig)
%     sim_T
    Delta = 0;
    Gamma = 0;
    Theta = 0;
    Vega = 0;
    Rho = 0;
    for month_i = 1:current_month
        rfr_T = T - (month_i - 1) * dT;
        r = RiskFreeRateInterpolation(rfr_T);
        
        LongAmount = LongCallAmounts(1, month_i);
        ShortAmount = ShortCallAmounts(1, month_i);
        LongStrike = LongCallStrikes(1, month_i);
        ShortStrike = ShortCallStrikes(1, month_i);
        
        LongValue =  LongAmount * bsm_call(r_stock, q, S0, LongStrike, sim_T, sig);
        ShortValue = ShortAmount * bsm_call(r_stock, q, S0, ShortStrike, sim_T, sig);
        
%         Can you just add these up?
        [delta, gamma, theta, vega, rho] = bsm_call_greeks(r_stock, q, S0, LongStrike, sim_T, sig);
        Delta = Delta + delta;
        Gamma = Gamma + gamma;
        Theta = Theta + theta;
        Vega = Vega + vega;
        Rho = Rho + rho;
        
%         Is the minus here correct?
        [delta, gamma, theta, vega, rho] = bsm_call_greeks(r_stock, q, S0, ShortStrike, sim_T, sig);
        Delta = Delta - delta;
        Gamma = Gamma - gamma;
        Theta = Theta - theta;
        Vega = Vega - vega;
        Rho = Rho - rho;
        
        
        BondValue = BondValue + BondFaceValues(1, month_i)*exp(-r * sim_T);
        OptionValue = OptionValue + (LongValue - ShortValue);
    end
    present_value = OptionValue + BondValue;
end


    




