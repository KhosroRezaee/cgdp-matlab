function [cal_idx, ts_idx] = split_val_for_calibration(val_idx, y, cfg)
%SPLIT_VAL_FOR_CALIBRATION Stratified split of validation indices into:
%  - ts_idx: for temperature scaling
%  - cal_idx: for conformal calibration
%
% Ensures each subset has at least cfg.model.min_pos_* positives when possible.

rng(cfg.model.seed + 123);
val_idx = val_idx(:);

pos = val_idx(y(val_idx)==1);
u   = val_idx(y(val_idx)==0);

pos = pos(randperm(numel(pos)));
u   = u(randperm(numel(u)));

% Target sizes ~50/50, but enforce minimum positives in each subset
npos = numel(pos);
minCal = min(cfg.model.min_pos_cal, npos);
minTS  = min(cfg.model.min_pos_ts,  max(0, npos-minCal));

% Default: half of positives to ts, rest to cal, then adjust to satisfy minima
npos_ts = floor(npos/2);
npos_ts = max(npos_ts, minTS);
npos_ts = min(npos_ts, npos - minCal);

npos_cal = npos - npos_ts;

ts_pos  = pos(1:npos_ts);
cal_pos = pos(npos_ts+1:end);

% Split unlabeled roughly 50/50 by count
nu = numel(u);
nu_ts = floor(nu/2);
ts_u = u(1:nu_ts);
cal_u = u(nu_ts+1:end);

ts_idx  = [ts_pos; ts_u];
cal_idx = [cal_pos; cal_u];

% Shuffle within each
ts_idx  = ts_idx(randperm(numel(ts_idx)));
cal_idx = cal_idx(randperm(numel(cal_idx)));

end
