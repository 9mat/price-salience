function [ theta, fval, exitFlag ] = fmincustom(routine, f, x0)
%FMINCUSTOM Summary of this function goes here
%   Detailed explanation goes here

if strcmp(routine, 'fminsearch')
    options = optimset;
    options = optimset(options, 'TolX', 1e-8);
    options = optimset(options, 'MaxIter', 100);
    options = optimset(options, 'Display', 'iter');
    
    [theta, fval, exitFlag] = fminsearch(f, x0, options);
    
elseif strcmp(routine, 'fminunc')
   options = optimoptions('fminunc');
   options.MaxIter = 10000;
   options.Display = 'iter';
   options.TolX = 1e-8;
   options.FinDiffType = 'central';
   options.MaxFunEvals = 10000;
   
   [theta, fval, exitFlag] = fminunc(f, x0, options);
   
elseif strcmp(routine, 'optitoolbox')
    lb = ones(size(x0)) *(-5);
    ub = ones(size(x0)) *5;
    lb(1) = -30; ub(1) = 30;
    
    options = optiset('maxiter', 1000);
    options = optiset(options, 'tolafun', 1e-6);
    options = optiset(options, 'tolrfun', 1e-6);
    options = optiset(options, 'display', 'iter');
    options = optiset(options, 'maxtime', 1e3);
    options = optiset(options, 'solver', 'ipopt');
    
    Opt = opti('fun', f, 'bounds', lb, ub, 'x0', x0, 'options', options);
    [theta, fval, exitFlag] = solve(Opt);
end

end

