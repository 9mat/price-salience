function f = nloglf(theta, Data, n)

likelihood = choiceProb(theta, Data, n);

if any(isnan(likelihood))
    f = 1e30;
    return;
end

loglf = log(likelihood);
loglf(likelihood<1e-100) = -1e6;
f = -sum(loglf);
end