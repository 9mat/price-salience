function f = nloglf(theta, Data, n)
likelihood = choiceProb(theta, Data, n);
loglf = log(likelihood);

if any(isnan(loglf))
    f = 1e20*(1+norm(theta)^2);
    return;
end

f = -sum(loglf);
end