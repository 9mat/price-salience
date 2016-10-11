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
    
    funcs.objective = f;
    funcs.gradient = @(x) finDiff(x, f);
    
    options.lb = lb;
    options.ub = ub;
    options.ipopt.hessian_approximation = 'limited-memory';
    
    % other options are in ipopt.opt
    
    [theta, info] = ipopt(x0, funcs, options);
    fval = info.eval.objective;
    exitFlag = info.status;
end

end

