function run_paper_simulation(pos_xlsx, neg_xlsx)
%RUN_PAPER_SIMULATION Backward-compatible entry point.
cfg = cfg_default();
run_paper_simulation_with_cfg(pos_xlsx, neg_xlsx, cfg);
end
