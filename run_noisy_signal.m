dataStruct = importdata('../data/data_new.csv');

data = dataStruct.data;
header = dataStruct.colheaders;

treat   = data(:,2)+1;
choice  = data(:,3);
price   = data(:,4:5);
demogr  = data(:,6:5);
X       = data(:,8:end);
X(:,end+1) = 1; % const

n.demogr = size(demogr, 2);
n.X      = size(X, 2);
n.choice = 3;
n.obs    = size(choice, 1);
n.treat  = max(treat);

n.gamma  = n.demogr;
n.beta   = n.X*(n.choice - 1);
n.sigma  = n.choice*(n.choice-1)/2;


price = price(:, 1:n.choice-1);

% matrix to change the base alternative (from the default 1 to i)
M = zeros(n.choice-1,n.choice-1,n.choice);
M(:,:,1) = eye(n.choice-1);
for i = 2:n.choice
    mm = eye(n.choice-1);
    mm(:, i-1) = -1;
    M(:,:,i) = mm;
end

n.M      = M;

Data.treat  = treat;
Data.choice = choice;
Data.price  = price;
Data.demogr = demogr;
Data.X      = X;

%% Estimation
theta0 = [
    -10; 
    % demogr
    %0; 0; 
  
    % beta ethanol
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;
    
    % beta midgrade gasoline
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;
    
    % sigma control
    0; 1; 

    % sigma treat1
    1; 0; 1
    
    % sigma treat2
    1; 0; 1
    ];

options1 = optimset('MaxIter', 2000, 'Display', 'iter');
options = optimoptions('fminunc', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'Display', 'iter', 'FinDiffType', 'central');
[theta] = fminsearch(@(x) nloglf(x,Data,n), theta0, options1);
[theta] = fminunc(@(x) nloglf(x,Data,n), theta, options);
[theta] = fminsearch(@(x) nloglf(x,Data,n), theta, options1);
[theta] = fminunc(@(x) nloglf(x,Data,n), theta, options);
[theta] = fminsearch(@(x) nloglf(x,Data,n), theta, options1);
[theta] = fminunc(@(x) nloglf(x,Data,n), theta, options);
[theta] = fminsearch(@(x) nloglf(x,Data,n), theta, options1);

%%

cov = mlecov(theta, Data, 'nloglf', @(x,d,a,b) nloglf(x,d,n));

%%
se = sqrt(diag(cov));
t = theta./se;

%%

rowheader = 'price';
for i = 1:n.demogr
    rowheader = [rowheader ' price*' header{5+i}];
end

for i=2:n.choice
    for j=1:n.X-1
        rowheader = [rowheader ' ' header{5+n.demogr+j} num2str(i)];
    end
    rowheader = [rowheader ' const' num2str(i)];
end

nn = 1 + n.gamma + n.beta;
printmat([theta(1:nn) se(1:nn) t(1:nn)], 'results', rowheader, 'est se t');

params = getParams(theta, n);

for i=1:n.treat
    fprintf('Variance matrix for treatment group %d\n', i);
    display(params.S(:,:,i)*params.S(:,:,i)');
end

%% bootstrap
% rng('default');
% nbstr = 1000;
% 
% theta_bstr = zeros(size(theta,1), nbstr);
% exit_flag = zeros(nbstr,1);
% options2 = optimset('MaxIter', 200, 'Display', 'iter');
% for i = 1:nbstr
%     bstr = randi(n.obs, n.obs, 1);
%     Data2.treat     = treat(bstr);
%     Data2.choice    = choice(bstr);
%     Data2.price     = price(bstr,:);
%     Data2.demogr    = demogr(bstr,:);
%     Data2.X         = X(bstr,:);
%     Data2.M         = M;
%     
%     theta0 = fminsearch(@(x) nloglf(x, Data2, n), theta, options2);
%     [theta_bstr(:,i),~,exit_flag(i)] = fminunc(@(x) nloglf(x, Data2, n), theta0, options);
% end

%% Average Marginal Effect

priceME = marginalPriceEffect(theta, Data, n);
varPriceME = covf(theta, @(x) marginalPriceEffect(x, Data, n), cov, size(priceME));

genderME = marginalXEffect(theta, Data, n, 1,1);
varGenderME = covf(theta, @(x) marginalXEffect(x, Data, n, 1, 1), cov, size(genderME));

ageME = zeros(n.choice, 3);
varAgeME = zeros(n.choice, 3);
for group = 1:3
    ageME(:,group) = marginalXEffect(theta, Data, n, group + 1, 2:4);
    varAgeME(:,group) = covf(theta, @(x) marginalXEffect(x, Data, n, group + 1, 2:4), cov, size(ageME(:,group)));
end

educationME = zeros(n.choice, 2);
varEducationME = zeros(n.choice, 2);
for group = 1:2
    educationME(:,group) = marginalXEffect(theta, Data, n, group + 4, 5:6);
    varEducationME(:,group) = covf(theta, @(x) marginalXEffect(x, Data, n, group + 4, 5:6), cov, size(varEducationME(:,group)));
end

cityME = zeros(n.choice, 4);
varCityME = zeros(n.choice, 4);
for group = 1:4
    cityME(:,group) = marginalXEffect(theta, Data, n, group + 6, 7:10);
    varCityME(:,group) = covf(theta, @(x) marginalXEffect(x, Data, n, group + 6, 7:10), cov, size(cityME(:,group)));
end

allME = [priceME genderME' ageME educationME cityME]';
allSeME = sqrt([varPriceME varGenderME' varAgeME varEducationME varCityME])';
allTME = allME./allSeME;

%% Counterfactual
probcf = cf(theta, n);
se_cf = sqrt(covf(theta, @(x) cf(x,n), cov, size(probcf)));
ci_l = probcf - 1.96*se_cf;
ci_u = probcf + 1.96*se_cf;

nprices = 1000;
price2 = linspace(0.206614, 0.4285888, nprices)'; % ethanol price 1%-99%-tile

fig=figure;
h=fill([price2; flipud(price2)] , [ci_l; flipud(ci_u)], 0.8*[1 1 1]); hold on;
set(h,'EdgeColor', 'None');
plot(price2, probcf, 'k');
xlabel('adjusted fuel price per km traveled (R$/km)');
ylabel('difference in prob choosing ethanol')
saveTightFigure(fig,'../report/cf.pdf');

%% Counterfactual 1.5 (comment out the last line in cf.m)
% prob = cf(theta,n);
% nprices = 1000;
% price2 = linspace(0.206614, 0.4285888, nprices)'; % ethanol price 1%-99%-tile
% fig1 = figure;
% plot(price2, prob(:,4), price2, prob(:,5));
% legend('Control group', 'Treatment group 1');
% xlabel('price of ethanol per km travelled (R$/km)');
% ylabel('prob of chossing ethanol');
% saveTightFigure(fig1,'../report/cfprob.pdf');

%% Counterfactual 2
fig2=figure;
prob1 = cf2(theta,n,1);
prob2 = cf2(theta,n,0.8) - prob1;
prob3 = cf2(theta,n,0.7) - prob1;
prob4 = cf2(theta,n,0.6) - prob1;
prob5 = cf2(theta,n,0.5) - prob1;
plot(price2,prob2(:,2), price2,prob3(:,2), price2,prob4(:,2), price2,prob5(:,2));
legend('0.8', '0.7', '0.6', '0.5');
xlabel('adjusted fuel price per km travelled (R$/km)');
ylabel('difference in prob choosing ethanol');
saveTightFigure(fig2,'../report/cf2.pdf');
