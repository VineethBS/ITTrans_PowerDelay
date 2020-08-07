%% Plot clambda as a function of lambda

H = [0.1, 1];
ph = [0.6, 0.4];

% H = [0.5, 1];
% ph = [0.5, 0.5];

% H = [0.75];
% ph = [1];

% H = [0.5, 1];
% ph = [0.9, 0.1];

% H = [1];
% ph = [1];

% H = [0.1, 0.5, 1];
% ph = [0.3, 0.4, 0.3];

% S = 0:0.1:2;
Smax = 3;

lambda = 0:0.0005:Smax;

clambda = zeros(1, length(lambda));

for i = 1:length(lambda);
    [clambda(i), ~] = minavgpowerreal(lambda(i), H, Smax, ph);
end

plot(lambda, clambda);
