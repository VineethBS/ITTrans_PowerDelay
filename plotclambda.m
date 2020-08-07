%% Plot clambda as a function of lambda

% H = [0.5, 1];
% ph = [1/3, 2/3];

% H = [0.5, 1];
% ph = [0.5, 0.5];

% H = [0.5, 1];
% ph = [0.7, 0.3];

% H = [sqrt(0.625)];
% ph = [1];

% H = [0.5, 1];
% ph = [0.9, 0.1];

% H = [1];
% ph = [1];

% H = [0.0001, 1, 2];
% ph = [0.2, 0.4, 0.4];

% H = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1];
% ph = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

% H = [0.1, 0.5, 1];
% ph = [0.3, 0.4, 0.3];

H = [0.1, 1];
ph = [0.6, 0.4];

% S = 0:0.1:2;
Smax = 3;
S = 0:Smax;

lambda = 0:0.005:Smax;
clambda = zeros(1, length(lambda));

for i = 1:length(lambda);
    [clambda(i), ~] = minavgpower(lambda(i), H, S, ph);
end

% subplot(2, 1, 1);
% plot(lambda, clambda, 'k');
% subplot(2, 1, 2);
% plot(lambda(1:end - 1), diff(clambda));

plot(lambda, clambda, 'k');