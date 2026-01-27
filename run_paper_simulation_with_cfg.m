function run_paper_simulation_with_cfg(pos_xlsx, neg_xlsx, cfg)
%RUN_PAPER_SIMULATION_WITH_CFG Paper-style end-to-end run on your local core dataset.
%
% This is a standalone function (not a subfunction), so it can be called from main.m.

if nargin < 3
    cfg = cfg_default();
end
if nargin < 2
    error('Usage: run_paper_simulation_with_cfg(pos_xlsx, neg_xlsx, cfg)');
end

if ~isfolder(cfg.out_dir); mkdir(cfg.out_dir); end

%% 1) Load + preprocess (no external sequences)
[seqs, y, meta] = load_core_dataset(pos_xlsx, neg_xlsx, cfg);

%% 2) Cluster split (<=40%% identity) and split clusters 80/10/10
splits = make_cluster_split(seqs, y, cfg);
save(fullfile(cfg.out_dir,'splits.mat'), 'splits', 'meta', 'cfg');

%% 3) Featurize
X = build_features(seqs, cfg);
tr = splits.train_idx; va = splits.val_idx; te = splits.test_idx;

%% 4) Train activity head with nnPU on Train
model = struct();
model.cfg = cfg;

rng(cfg.model.seed);
model.activity = train_nnpu_logreg(X(tr,:), y(tr), cfg.model.pi, cfg);

%% 5) Temperature scaling + 6) Conformal gate
[cal_idx, ts_idx] = split_val_for_calibration(va, y, cfg);

logits_ts = predict_logit(model.activity, X(ts_idx,:));
T = fit_temperature(logits_ts, y(ts_idx)); % proxy labels: P=1, U=0
model.temperature = T;

p_cal = sigmoid(predict_logit(model.activity, X(cal_idx,:)) / T);
gate = fit_conformal_gate(p_cal, y(cal_idx), cfg.model.alpha, cfg.model.min_p_act);
model.conformal = gate;

save(fullfile(cfg.out_dir,'model.mat'), 'model');

%% 7) Evaluate on test (proxy evaluation: P vs U)
p_te = sigmoid(predict_logit(model.activity, X(te,:)) / model.temperature);
metrics = evaluate_binary_metrics(y(te), p_te);

json_path = fullfile(cfg.out_dir, 'test_metrics.json');
fid=fopen(json_path,'w'); fprintf(fid, '%s', jsonencode(metrics)); fclose(fid);

fprintf('Test AUROC=%.3f  AUPRC=%.3f  (proxy labels P vs U)\n', metrics.auroc, metrics.auprc);
fprintf('Saved metrics to %s\n', json_path);

%% 8) Run CTCM‑Neo generation (paper hyperparams), using Train positives as seeds
train_pos = tr(y(tr)==1);
seed_seqs = seqs(train_pos);

for r = 1:cfg.ctcm.seeded_runs
    rng(100+r);
    archive = ctcm_neo_generate(seed_seqs, model, cfg);
    out_csv = fullfile(cfg.out_dir, sprintf('ctcm_archive_run%d.csv', r));
    write_archive_csv(out_csv, archive);
    fprintf('Run %d: wrote %s (%d candidates)\n', r, out_csv, numel(archive));
end

end
