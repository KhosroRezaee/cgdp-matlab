function gate = fit_conformal_gate(p_cal, y_cal, alpha, min_p_act)
%FIT_CONFORMAL_GATE Split conformal threshold on positive calibration points.
% Nonconformity score: s = 1 - p.
% Threshold q = (1-alpha)-quantile of scores among positives.
% Accept if p >= max(min_p_act, 1-q).

pos = (y_cal==1);
if sum(pos) < 3
    warning('Too few positives in calibration; using all points as fallback.');
    pos = true(size(y_cal));
end

scores = 1 - p_cal(pos);

q = quantile_higher(scores, 1 - alpha);
t_conf = 1 - q;

gate = struct();
gate.alpha = alpha;
gate.q = q;
gate.t_conf = t_conf;
gate.t_act = max(min_p_act, t_conf);
end

function q = quantile_higher(x, qlevel)
% Conservative "higher" quantile like conformal uses.
x = sort(x(:));
n = numel(x);
k = ceil(qlevel * (n+1));
k = min(max(k,1), n);
q = x(k);
end
