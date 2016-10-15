function [ grad ] = finDiff( theta, f )
%FINDIFF Summary of this function goes here
%   Detailed explanation goes here

f0 = f(theta);
grad = zeros(numel(f0), numel(theta));

e = 1e-5;

for i=1:numel(theta)
    theta_l = theta;
    theta_u = theta;
    
    if abs(theta(i)) > e*1e-2
        theta_l(i) = theta(i)*(1-e);
        theta_u(i) = theta(i)*(1+e);
    else
        theta_l(i) = theta(i)-e;
        theta_u(i) = theta(i)+e;
    end
    
    f_l = f(theta_l);
    f_u = f(theta_u);
    
    grad(:, i) = (f_u(:) - f_l(:))/(theta_u(i) - theta_l(i));
    grad(isnan(grad(:,i)),i) = 1e10*(-theta(i));
end

end

