function [ prob ] = cf( theta, n)
%CF Summary of this function goes here
%   Detailed explanation goes here

nprices = 1000;

sample_demogr = [
    0; % expensive car
    0; % heavy usage
    ];

sample_X = [
    0; % male
    
    0; % not top 75 pctile car price
    0; % not top 75 pctile car usage
    
    1; % age 25-40;
    0; % age 40-65
    0; % age > 65
    
    1; % some secondary
    0; % some college
    
    1; % SP
    0; % CTB
    0; % BH
    0; % REC
    
    1 % const
    ];

price1 = 0.2701633*ones(nprices,1); % regular gasoline median price
price3 = 0.2835509*ones(nprices,1); % midgrade gasoline median price
price2 = linspace(0.206614, 0.4285888, nprices)'; % ethanol price 1%-99%-tile
price = [price2 - price1, price3 - price1];

% order: choice -> treat -> price
Data.choice     = reshape(repmat(1:n.choice,n.treat*nprices, 1), n.choice*n.treat*nprices, 1);
Data.treat      = repmat(reshape(repmat(1:n.treat,nprices,1), nprices*n.treat, 1), n.choice, 1);
Data.price      = repmat(price, n.choice*n.treat, 1);
Data.demogr     = repmat(sample_demogr', n.choice*n.treat*nprices, 1);
Data.X          = repmat(sample_X', n.choice*n.treat*nprices, 1);

prob = choiceProb(theta, Data, n);
prob = reshape(prob, nprices, n.treat*n.choice);
prob = prob(:,5) - prob(:,4);
end

