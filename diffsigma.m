function [ ds ] = diffsigma( theta, n )
%DIFFSIGMA Summary of this function goes here
%   Detailed explanation goes here

params = getParams(theta, n);

for i = 1:n.treat
    v{i} = params.S(:,:,i)*params.S(:,:,i)';
end

ds = zeros(n.treat-1);

for i = 1:n.treat-1
    ds(:,i) = diag(v{i+1}-v{1});
end

ds = ds(:);
end

