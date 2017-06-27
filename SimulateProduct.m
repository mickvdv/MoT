% product specifications
premium = 150;
global payments_per_year;
payments_per_year = 12;

N = 1000;

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
sig = 0.219;
q = 0.0203;

% pricing
global r_stock;
r_stock = q + 0.062;
month_duration = 21; % trading days per month

PayOffs = zeros(N,1);

for sim_i = 1:N
    StockPath = GenerateStockPath(S0, T, (1/(month_duration*payments_per_year)), (r_stock-q), sig);

    LongCallStrikes = zeros(1, T/dT);
    LongCallAmounts = zeros(1, T/dT);
    ShortCallAmounts = zeros(1, T/dT);
    ShortCallStrikes = zeros(1, T/dT);
%     Deltas = zeros(1, T/dT);
%     Gammas = zeros(1, T/dT);
%     Thetas = zeros(1, T/dT);
%     Vegas = zeros(1, T/dT);
%     Rhos = zeros(1, T/dT);

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

    S0 = StockPath(1, 252*10);

    [present_value, Delta, Gamma, Theta, Vega, Rho] = CalculateProductPrice(S0, LongCallStrikes, LongCallAmounts, ShortCallStrikes, ShortCallAmounts, BondFaceValues, T - (day_i/(month_duration*payments_per_year)));

    PayOffs(sim_i) = present_value / sum(BondFaceValues);
end

mean(PayOffs)


% figure
% subplot(4,1,1)
% plot(StockPath);
% yyaxis right
% xlabel('Time (in days)')
% plot(values);
% legend('Stock path', 'Product value');
% 
% subplot(4,1,2)
% hold on
% plot(Deltas)
% plot(Gammas)
% hold off
% xlabel('Time (in days)')
% legend('Delta', 'Gamma');
% 
% subplot(4,1,3)
% hold on
% plot(Vegas)
% plot(Rhos)
% hold off
% xlabel('Time (in days)')
% legend('Vega', 'Rho');
% 
% subplot(4,1,4)
% plot(Thetas)
% xlabel('Time (in days)')
% legend('Theta');

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
        
% %         Can you just add these up?
%         [delta, gamma, theta, vega, rho] = bsm_call_greeks(r_stock, q, S0, LongStrike, sim_T, sig);
%         Delta = Delta + delta;
%         Gamma = Gamma + gamma;
%         Theta = Theta + theta;
%         Vega = Vega + vega;
%         Rho = Rho + rho;
%         
% %         Is the minus here correct?
%         [delta, gamma, theta, vega, rho] = bsm_call_greeks(r_stock, q, S0, ShortStrike, sim_T, sig);
%         Delta = Delta - delta;
%         Gamma = Gamma - gamma;
%         Theta = Theta - theta;
%         Vega = Vega - vega;
%         Rho = Rho - rho;
        
        
        BondValue = BondValue + BondFaceValues(1, month_i)*exp(-r * sim_T);
        OptionValue = OptionValue + (LongValue - ShortValue);
    end
    present_value = OptionValue + BondValue;
end


    




