format long g

% simulation parameters
N = 500; % number of paths in Monte Carlo


% for loading the basket
filename = 'basket.xlsx';

% Load the most recent risk free interpolation
global risk_free_rate_interpolation;
load('rfr.mat', 'risk_free_rate_interpolation');

r_stock = RiskFreeRateInterpolation(3/12);

% product specifications
premium = 150;
global payments_per_year;
payments_per_year = 12;
T = 10;
dT = 1/payments_per_year;

% basket parameters
S0 = 68.63;
sig = 0.1567;
q = 0.0203;

% CapRates = 1:0.025:1.5;
CapRates = [1.10];
% ParticipationRates = 0.2:0.025:1.2;
ParticipationRates = [0.95];

ProtectionRates = zeros(size(CapRates,2), size(ParticipationRates, 2));

for cap_i = 1:size(CapRates,2)
    cap_rate = CapRates(1, cap_i);
    for participation_i = 1:size(ParticipationRates, 2)
        participation_rate = ParticipationRates(participation_i);
        % calculate the call prices
        CallPrices = zeros(1, T*payments_per_year);
        CallAmounts = zeros(1, T*payments_per_year);


        for month_i = 1:(T/dT)
            sim_T = T - (month_i - 1) * dT;
%           sim_T = time left till maturity
            r = RiskFreeRateInterpolation(sim_T);


    %         get the expected values of the calls
            ExpectedLongCallValue = bsm_call(r, q, S0, S0, sim_T, sig);
            K_short = S0 * cap_rate ^ (sim_T);
            ExpectedShortCallValue = bsm_call(r, q, S0, K_short, sim_T, sig);

            amount_of_stock = (premium * participation_rate) / S0;


    %         save prices and amounts
            CallPrices(1, month_i) = ExpectedLongCallValue * amount_of_stock - ExpectedShortCallValue * amount_of_stock;
            CallAmounts(1, month_i) = amount_of_stock;

    %         calculated the expected value for next month
            StockPaths = SimulateStockPaths(S0, dT, dT, r_stock - q, sig, N);
            S0 = ExpectedValueFromStockPaths(StockPaths, dT, dT, N);     
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
        protection_rate = Protection / (premium * payments_per_year * T);
        ProtectionRates(cap_i, participation_i) = protection_rate;
    end
end

save('bonds.mat', 'BondFaceValues');

hold on
surf(ParticipationRates, CapRates, ProtectionRates);
title('Protection Rate');
ylabel('Cap Rate');
xlabel('Participation Rate');
hold off



