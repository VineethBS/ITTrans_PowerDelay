%% Policyiteration

function [aq, ac, pol_new] = politer(Beta, H, pH, S, Amax, lambda)
lH = length(H);
lS = length(S);
Smax = S(lS);

if Beta < 1000
    Qmax = ceil(Beta + 1000);
else
    Qmax = ceil(2 * Beta);
end

% Qmax = ceil(Beta);

% Qmax = ceil(5 * Beta);
pA = binopdf(0:Amax, Amax, lambda/Amax);

disp(datestr(now));
% the inital policy vector
pol_old = Smax * ones(Qmax, lH);
% a simple estimate
for q = 1:Qmax
    pol_old(q, :) = min([floor(q/Beta), Smax]);
end
pol_old(1, :) = 0;
pol_new = pol_old; %% NOTE : Stores batch size, not the index !!!

pos = @(q, h) ((q - 1) * lH + (h - 1) + 1); % finds the position in the single dim array
ql = @(q) (q - 1); % queue length corresponding to the index

maxiterations = 20;

for iter = 1:maxiterations
    % fprintf(1,'Iter %u\n', iter);
    % Setting up the sparse coefficient matrix
    rows = [];
    cols = [];
    vals = [];
    for q = 1:Qmax
        for h = 1:lH
            rows = [rows, pos(q,h)];
            cols = [cols, pos(q,h)];
            vals = [vals, 1];
            for nh = 1:lH
                for a = 0:Amax
                    n = min([max([ql(q) - pol_old(q, h), 0]) + a, Qmax - 1]);
                    rows = [rows, pos(q,h)];
                    cols = [cols, pos(n + 1,nh)];
                    vals = [vals, -pH(nh) * pA(a + 1)];
                end
            end
        end
    end
    
    coeff = sparse(rows, cols, vals);
    coeff(:,1) = 1;
    
    costvector = zeros(Qmax * lH, 1);
    
    for q = 1:Qmax
        for h = 1:lH
            costvector(pos(q,h)) = ql(q) + Beta * powercost( H(h), pol_old(q, h) );
        end
    end
    
    J = coeff\costvector;
    g = J(1);
    J(1) = 0;

    % Do policy improvement
    for q = 1:Qmax
        for h = 1:lH
            max_sforq = min(ql(q), Smax);
            costs = zeros(1, max_sforq + 1);
            for s = 0:max_sforq
                costs(s + 1) = ql(q) + Beta * powercost(H(h), s) - g;
                for nh = 1:lH
                    for a = 0:Amax
                        n = min(ql(q) - s + a, Qmax - 1);
                        costs(s + 1) = costs(s + 1) + pH(nh) * pA(a + 1) * J(pos(n + 1, nh));
                    end
                end
            end
            [~, pol_new(q, h)] = min(costs);
        end
    end
    pol_new = pol_new - 1;
    % Check if policies have converged
    if max(max(abs(pol_old - pol_new))) == 0
        break;
    end
    pol_old = pol_new;     
end

% Evaluate the good policy but only the queue length
% Setting up the sparse coefficient matrix
rows = [];
cols = [];
vals = [];

for q = 1:Qmax
    for h = 1:lH
        rows = [rows, pos(q,h)];
        cols = [cols, pos(q,h)];
        vals = [vals, 1];
        for nh = 1:lH
            for a = 0:Amax
                n = min([max([ql(q) - pol_old(q, h), 0]) + a, Qmax - 1]);
                rows = [rows, pos(q,h)];
                cols = [cols, pos(n + 1,nh)];
                vals = [vals, -pH(nh) * pA(a + 1)];
            end
        end
    end
end

coeff = sparse(rows, cols, vals);
for q = 1:(Qmax * lH)
    coeff(q, 1) = 1;
end

costvector = zeros(1, Qmax * lH);

for q = 1:Qmax
    for h = 1:lH
        costvector(pos(q,h)) = ql(q);
    end
end

costvector = costvector';
J = coeff\costvector;
aq = J(1);
J(1) = 0;

ac = (g - aq)/Beta;

disp(datestr(now))
