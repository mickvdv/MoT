
payments_per_year = 12;

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

figure
plot(x,risk_free_rate,'o')
hold on
plot(x,risk_free_rate_interpolation)
xlabel('T (in months)');
ylabel('Risk Free Rate (in %)');
hold off
risk_free_rate_interpolation = risk_free_rate_interpolation / 100;
save('rfr.mat', 'risk_free_rate_interpolation');

