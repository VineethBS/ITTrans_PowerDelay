%% Example for minimum average power

function [clambda, ps] = minavgpower(lambda, H, S, ph)

ps0 = (1/length(S)) * ones(length(H), length(S));

[ps, clambda] = fmincon(@avgpower, ps0, [], [], [], [], zeros(length(H), length(S)), ones(length(H), length(S)), @mycon);

    function p = avgpower(ps)
        p = 0;
        for i = 1:length(H)
            p = p + ph(i) * sum(ps(i,:) .* (1.28 * (10.^(50 * S/200) - 1)/(H(i)^2)));
            % p = p + ph(i) * sum(ps(i,:) .* (2 * S/(H(i))));
        end
    end
    function [c, ceq] = mycon(ps)
        ceq = [ratecons(ps); pscons(ps)'];
        c = [];
    end
    function c = ratecons(ps)
        c = 0;
        for i = 1:length(H)
            c = c + ph(i) * sum(ps(i,:) .* S);
        end
        c = c - lambda;
    end

    function c = pscons(ps)
        c = [];
        for i = 1:length(H)
            c = [c, sum(ps(i,:)) - 1];
        end
    end
end
    
    
    
