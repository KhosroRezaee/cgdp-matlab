function cfg = cfg_default()
% Default configuration approximating the paper draft.
% Keep all knobs here so you can reproduce experiments.

%% Data / preprocessing
cfg.data.min_len = 8;
cfg.data.max_len = 30;
cfg.data.allow_ambiguous = false;   % only canonical 20 AA
cfg.data.dedup = true;

% Handling of sequences longer than max_len in provided files:
%  - "drop"        : remove them (paper-faithful if positives already in 8–30)
%  - "truncate"    : take first max_len residues
%  - "best_window" : choose the length=max_len window maximizing (charge - 0.5*|gravy| - 0.5*boman)
cfg.data.long_strategy = "best_window";

%% Physicochemical windows used for generation (paper)
cfg.windows.charge_min = +3;
cfg.windows.charge_max = +7;
cfg.windows.gravy_min  = -1.5;
cfg.windows.gravy_max  = +0.5;
cfg.windows.boman_max  = 1.5;

%% Clustering / split
cfg.cluster.use_cdhit = true;
cfg.cluster.cdhit_identity = 0.40; % <=40% identity
cfg.cluster.train_frac = 0.80;
cfg.cluster.val_frac   = 0.10;
cfg.cluster.test_frac  = 0.10;
cfg.cluster.seed = 7;

% Ensure enough positives land in val/test for calibration (especially when not using CD-HIT)
cfg.cluster.min_pos_val = 5;
cfg.cluster.min_pos_test = 5;

% CD-HIT executable (edit if needed)
cfg.cluster.cdhit_bin = "cd-hit";

%% Model: ConformaX‑PEP (MATLAB-native surrogate)
cfg.model.pi = 0.22;        % fixed PU prior (paper)
cfg.model.alpha = 0.10;     % conformal risk level
cfg.model.min_p_act = 0.80; % operating point threshold (paper ~0.78–0.80)
cfg.model.use_hemo = false; % set true if you add hemolysis labels later
cfg.model.max_p_hemo = 0.20;

% Feature choices
cfg.model.embedding = "dipeptide"; % "dipeptide" or "none" or "precomputed"
cfg.model.zscore = true;

% nnPU logistic regression training
cfg.model.lr = 0.05;
cfg.model.epochs = 2500;
cfg.model.batch_size = 64;
cfg.model.weight_decay = 1e-4;
cfg.model.nnpu_beta = 0.0;     % optional non-negativity slack (0 = strict)
cfg.model.seed = 11;

% Temperature scaling
cfg.model.temp_init = 1.0;

% Val split controls (small datasets need stratification)
cfg.model.min_pos_cal = 3;  % minimum positives in conformal calibration subset
cfg.model.min_pos_ts  = 2;  % minimum positives in temperature-scaling subset

%% CTCM‑Neo (paper-like hyperparameters)
cfg.ctcm.population = 128;
cfg.ctcm.num_tribes = 8;                % 8 tribes x 16 members = 128
cfg.ctcm.archive = 16;
cfg.ctcm.seeded_runs = 5;
cfg.ctcm.proposals_per_run = 50000;

cfg.ctcm.p_sub = 0.75;
cfg.ctcm.p_ins = 0.125;
cfg.ctcm.p_del = 0.125;

% annealing / schedules
cfg.ctcm.T0 = 1.0;
cfg.ctcm.Tmin = 0.05;
cfg.ctcm.anneal_horizon = cfg.ctcm.proposals_per_run;
cfg.ctcm.burnin_frac = 0.25;  % coarse-to-fine switch point
cfg.ctcm.metropolis = true;

% diversity (simple)
cfg.ctcm.min_hamming_identity = 0.85; % prune near duplicates in archive

% novelty distance weight
cfg.ctcm.lambda_novel = 1.0;
cfg.ctcm.lambda_len = 0.05;
cfg.ctcm.lambda_penalty = 10.0;

%% Output
cfg.out_dir = "out";

end
