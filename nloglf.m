function f = nloglf(theta, Data, n)

likelihood = choiceProb(theta, Data, n);

if any(isnan(likelihood))
    f = 1e30;
    return;
end

likelihood(likelihood<1e-100) = 1e-100;
f = -sum(log(likelihood));
end