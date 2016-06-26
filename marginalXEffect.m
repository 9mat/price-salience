function [ ame ] = marginalXEffect( theta, Data, n, i, ii )
%MARGINALEFFECTS Summary of this function goes here
%   Detailed explanation goes here

% detect binary variable
values = unique(Data.X(:,i));

isbinary = (numel(values) == 2) & (max(values) == 1) & (min(values) == 0);

prob0 = zeros(n.obs, n.choice);
prob1 = zeros(n.obs, n.choice);

if isbinary
    Data.X(:, ii) = 0;
else
    d = Data.X(:, i)*1e-3;
    d(abs(d)<1e-6) = 1e-6;
    Data.X(:,i) = Data.X(:,i) - d;
end

for j = 1:n.choice
    Data.choice(:) = j;
    prob0(:,j) = choiceProb(theta, Data, n);
end
    
if isbinary
    Data.X(:,i) = 1;
else
    Data.X(:,i) = Data.X(:,i) + 2*d;
end    

for j = 1:n.choice
    Data.choice(:) = j;
    prob1(:,j) = choiceProb(theta, Data, n);
end

ame = prob1 - prob0;

if ~isbinary
    ame = bsxfun(@rdivide, ame, 2*d);
end

ame = squeeze(mean(ame, 1));

end

