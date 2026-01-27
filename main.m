function main(pos_xlsx, neg_xlsx)
%MAIN One-stop entry point for the paper-style simulation.
%
% Usage:
%   main('peptides.xlsx','negatives_200_generated.xlsx')

if nargin < 2
    error('Usage: main(pos_xlsx, neg_xlsx)');
end

% Add repo to MATLAB path
addpath(genpath(pwd));

cfg = cfg_default();

% Auto-disable CD-HIT if missing (Windows users commonly hit this)
if cfg.cluster.use_cdhit && ~command_exists(cfg.cluster.cdhit_bin)
    fprintf('[main] cd-hit not found; switching to MATLAB-only clustering fallback.\n');
    cfg.cluster.use_cdhit = false;
end

run_paper_simulation_with_cfg(pos_xlsx, neg_xlsx, cfg);
end
