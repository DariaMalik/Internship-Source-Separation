function y = kurtosis(x)

y = mean(x.^4) - (3 * mean(x.^2) * mean(x.^2));