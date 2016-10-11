function f = nloglf(theta, Data, n)
likelihood = choiceProb(theta, Data, n);

if any(isnan(likelihood))
    f = 1e20*(1+norm(theta)^2);
    return;
end

loglf = log(likelihood);
f = -sum(loglf);
end