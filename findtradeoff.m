%% Find tradeoff

% betas = [10, 100, 1000];
betas = [4000];

aq = zeros(1, length(betas));
ac = zeros(1, length(betas));

H = [0.1, 1];
pH = [0.6, 0.4];

lambda = 0.80;
tradeoff_file = 'Tradeoff_H(0.1,1)_pH(0.6,0.4)_lambda_0.8.mat';

% lambda = 2.395;
% tradeoff_file = 'Tradeoff_H(0.1,0.5,1)_pH(0.3,0.4,0.3)_lambda_2.395.mat';
% 
% lambda = 2.405;
% tradeoff_file = 'Tradeoff_H(0.1,0.5,1)_pH(0.3,0.4,0.3)_lambda_2.405.mat';

S = [0,1,2,3];
Amax = 5;

for i = 1:length(betas)
    Beta = betas(i);
    disp(sprintf('\nBeta %f', Beta));
    [aq(i), ac(i), ~] = politer(Beta, H, pH, S, Amax, lambda);
    % [aq(i), ac(i), p] = politer(Beta, H, pH, S, Amax, lambda);
    save(tradeoff_file);
    % save policy if needed
    % fname = sprintf('%u_policy.mat',Beta);
    % save(fname, 'p');
end