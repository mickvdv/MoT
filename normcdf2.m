%

function y = normcdf2(x, mu, sigma)

f = @(u,o,x) 0.5*(1+erf((x-u)/sqrt(2*o^2)));

y = feval(f, mu, sigma, x);


