function T = fit_temperature(logits, y)
%FIT_TEMPERATURE Temperature scaling (single scalar) by minimizing Brier loss.
% y: proxy labels (1 for P, 0 for U)

% optimize log(T) so T>0
obj = @(logT) brier(sigmoid(logits ./ exp(logT)), y);
logT0 = log(1.0);
opts = optimset('Display','off');
logT = fminsearch(obj, logT0, opts);
T = exp(logT);
end

function b = brier(p, y)
b = mean((p - y).^2);
end
