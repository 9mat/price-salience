function f = nloglf(theta, Data, n)
likelihood = choiceProb(theta, Data, n);

loglf = log(likelihood);
loglf(likelihood<1e-100) = 1e5*loglf(likelihood<1e-100);
loglf(isnan(loglf)) = -1e6*(1+rand);
f = -sum(loglf);
end