# CGDP / CTCM‑Neo + ConformaX‑PEP (MATLAB) — paper-style simulation on your core dataset (no external sequences)

This repo reproduces the *simulation pipeline* described in the paper draft you shared:
- Core data: 52 positives (P) + 200 unlabeled/background (U). No external sequence data; no positive-like; no external test.
- Preprocess: canonical AA, length 8–30, deduplication; compute descriptors (charge@pH7, GRAVY, Boman, length).
- Cluster split: CD-HIT at <=40% identity; split clusters into Train/Val/Test = 80/10/10.
- ConformaX‑PEP (simplified, MATLAB-native):
  - Feature encoder: PhysChem + optional embeddings (default: AA composition + dipeptide composition).
  - Activity head: nnPU logistic regression with fixed prior pi=0.22.
  - Temperature scaling on Val.
  - Split conformal gate (alpha=0.10) using calibration subset from Val.
  - Hemolysis head is **optional**: if you provide labels, it trains; otherwise gating skips hemolysis.
- CTCM‑Neo generator (paper hyperparams): population=128, archive=16, p_sub=0.75, p_ins=0.125, p_del=0.125, proposals=50k/run, seeded_runs=5
  - Variation ops: substitution/insertion/deletion with guided repair into windows.
  - Coarse-to-fine objective: novelty/feasibility → add calibrated p_active.
  - Accept/reject: Metropolis + conformal gate at end.

## Files
- `run_paper_simulation.m` : main entry point
- `cfg_default.m` : all knobs (windows, alpha, pi, CTCM hyperparams)
- `data/load_core_dataset.m`
- `featurize/*` : descriptors + feature vectors
- `cluster/*` : FASTA export + CD-HIT parsing + cluster split
- `model/*` : nnPU logistic reg, temperature scaling, conformal gate
- `ctcm/*` : CTCM‑Neo generator
- `utils/*` : helper functions

## Requirements
- MATLAB (R2020b+ recommended). No toolboxes strictly required.
- CD-HIT installed and on PATH **if you want exact clustering like paper**.
  - If CD-HIT is not available, the code can fall back to a naive hash-based “cluster” (NOT equivalent).
  - Paper-faithful: install CD-HIT and keep `cfg.cluster.use_cdhit = true`.

## Run
Open MATLAB at the repo folder and execute:
```matlab
main('peptides.xlsx','negatives_200_generated.xlsx');
% (or) run_paper_simulation('peptides.xlsx','negatives_200_generated.xlsx');
```
Outputs:
- `out/splits.mat` (indices, clusters)
- `out/model.mat` (weights, temperature, conformal threshold)
- `out/ctcm_archive_run*.csv` (generated peptides + scores)
- `out/test_metrics.json` (AUROC/AUPRC etc. on proxy labels P vs U)

## Notes about “paper fidelity”
- The paper describes a frozen ProtT5 backbone. This repo keeps the interface *pluggable*:
  - Default embeddings: AA composition + dipeptide composition (fully local, no model download).
  - If you want ProtT5 embeddings: replace `featurize/seq_embedding.m` to load your precomputed embeddings.
- Hemolysis head requires hemolysis labels (not present in your current core set). If you later provide a hemolysis-labeled table, set `cfg.model.use_hemo=true`.



## Windows note (your error)
If you see "'cd-hit' is not recognized..." on Windows, run `main(...)` which auto-falls back to a MATLAB-only clustering method. For exact paper clustering, install CD-HIT and ensure `cd-hit` is on PATH.

## Positive count note
If your positives file contains sequences longer than 30, the paper's strict filter would drop them. To keep your provided positive count while enforcing length<=30, this repo uses `cfg.data.long_strategy = "best_window"` by default (crops long sequences to the best 30-aa window by a simple physchem heuristic). Set it to "drop" for strict paper behavior.


### Fix in v3
If you previously saw `containers.Map` key/value mismatch, this is fixed by using CHAR keys in `featurize/aa_tables.m`.


### Entry points
- `main(...)` (recommended)
- `run_paper_simulation(...)` (wrapper)
- `run_paper_simulation_with_cfg(...)` (advanced: pass custom cfg)


### Fix in v5
MATLAB-only clustering fallback produced variable-shaped cluster index vectors; `make_cluster_split` now uses a robust flattener instead of `cell2mat`.


### Fix in v6
- Fixed random amino-acid sampling in CTCM (`random_aa` and swap helpers) to avoid `extractBetween` start>end.
- Added `cfg.cluster.min_pos_val/min_pos_test` and enforcement to prevent Val/Test from having too few positives (helps conformal calibration).


### Fix in v7
- Fixed MATLAB syntax error in CTCM (`'KR'(idx)` is invalid). Now uses a variable `opts='KR'` before indexing.
- Validation split for temperature vs conformal is now stratified to reduce 'Too few positives' warnings.

## Reference
This code accompanies the paper:
**Intelligent in silico prioritization of antimalarial peptide candidates under explicit physicochemical windows via de novo CTCM-Neo generation and conformal-gated calibrated classification**
Muhammad Aamir, Khosro Rezaee* and Maryam Saberi Anari

## Reference
This code accompanies the paper:
**Intelligent in silico prioritization of antimalarial peptide candidates under explicit physicochemical windows via de novo CTCM-Neo generation and conformal-gated calibrated classification**
Muhammad Aamir, Khosro Rezaee* and Maryam Saberi Anari
