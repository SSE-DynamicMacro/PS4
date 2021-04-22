function [foc] = transition_path(guess,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c = guess(:,1);
l = guess(:,2);
k = guess(:,3);
lambda = guess(:,4);

alpha = params(1);
beta = params(2);
delta = params(3);
theta = params(4);
A = params(5);
k_ss = params(6);

l_lead = [l(2:end);l(end)];
lambda_lead = [lambda(2:end);lambda(end)*beta];
beta = cumprod(repmat(beta,1,length(guess))');

kk = [k_ss;k(1:end-1)];

foc = NaN(length(guess),1);
foc(:,1) = lambda - beta .* alpha ./ c;
foc(:,2) = (beta .* (1 - alpha)) ./ (1 - l) - (lambda * (1 - theta) * A ^ (1 - theta)) .* kk .^ theta .* l .^ (-theta);
foc(:,3) = lambda_lead .* (theta * A ^ (1 - theta) * k .^ (theta - 1) .* l_lead .^ (1 - theta) + (1 - delta)) - lambda;
foc(:,4) = A ^(1 - theta) * kk .^ theta .* l .^ (1 - theta) + (1 - delta) * kk - c - k;
end

