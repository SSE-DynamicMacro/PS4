cd ../PS4
clear
clc

%%
N = 100;

A0 = 1; A1 = 1.1;

k_bar = 3; % k_bar = k/y
theta = 1/3;
delta = 0.1;

l_ss = 1/3;
y_ss = k_bar ^(theta/(1-theta)) * l_ss;
k_ss = k_bar * y_ss;
c_ss = y_ss - delta * k_ss;

r = theta * k_bar^(-1) - delta;
beta = 1/(1 + r);

alpha = (1-delta*k_bar) / ((1-theta) / l_ss - (1-theta) + (1-delta*k_bar));

lambda = alpha / c_ss;

%%
guess = [repmat([c_ss l_ss k_ss],N,1) cumprod([lambda repmat(beta,1,N-1)])'];
params = [alpha beta delta theta A1 k_ss];
[X] = fsolve(@transition_path,guess,optimset('Display','none'),params);

c = X(:,1);
p.consumption = [c_ss;c];
l  = X(:,2);
p.labor = [l_ss;l];
k = X(:,3);
p.capital = [k_ss;k];

y = theta * A1 ^ (1 - theta) * k .^ (theta - 1) .* l .^ (1 - theta);
p.output = [y_ss;y];

%%
eps=1e-5;

figure(1)
subplot(2,2,1)
plot(p.output,'LineWidth',1.5)
title('Output')
xlim([0 N])
ylim([min(p.output)*(1-eps) max(p.output)*(1+eps)])
subplot(2,2,2)
plot(p.consumption,'LineWidth',1.5)
title('Consumption')
xlim([0 N])
ylim([min(p.consumption)*(1-eps) max(p.consumption)*(1+eps)])
subplot(2,2,3)
plot(p.labor,'LineWidth',1.5)
title('Labor')
xlim([0 N])
ylim([min(p.labor)*(1-eps) max(p.labor)*(1+eps)])
subplot(2,2,4)
plot(p.capital,'LineWidth',1.5)
title('Capital')
xlim([0 N])
ylim([min(p.capital)*(1-eps) max(p.capital)*(1+eps)])