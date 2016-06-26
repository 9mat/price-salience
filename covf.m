function [ v ] = covf( theta, f, cov, sz )
%COVF Summary of this function goes here
%   Detailed explanation goes here

df_theta = finDiff(theta, f);
v = sum(df_theta.*(df_theta*cov'),2);
v = reshape(v, sz);

end

