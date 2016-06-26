function [ params ] = getParams( theta, n )
%GETPARAMS Summary of this function goes here
%   Detailed explanation goes here

offset = 0;

params.alpha = theta(1);
offset = offset + 1;

params.gamma = theta(offset + 1: offset + n.gamma);
offset = offset + n.gamma;

params.beta = theta(offset + 1: offset + n.beta);
params.beta = reshape(params.beta, n.X, n.beta/n.X);
offset = offset + n.beta;

params.S = zeros(n.choice - 1, n.choice - 1, n.treat);

% mark the lower triangular matrix
tril_mark = tril(ones(n.choice-1)) == 1;

mat_S = zeros(n.choice-1);

vec_S = theta(offset + 1: offset + n.sigma - 1);
mat_S(tril_mark) = [1; vec_S];
params.S(:,:,1) = mat_S;
offset = offset + n.sigma - 1;

for t = 2:n.treat
    vec_S = theta(offset + 1: offset + n.sigma);
    mat_S(tril_mark) = vec_S;
    params.S(:,:,t) = mat_S;
    offset = offset + n.sigma;
end

end

