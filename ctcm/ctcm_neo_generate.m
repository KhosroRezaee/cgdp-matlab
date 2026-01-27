function archive = ctcm_neo_generate(seed_seqs, model, cfg)
%CTCM_NEO_GENERATE Paper-like CTCM‑Neo discrete swarm for peptides.
% Returns an array of structs with fields:
%   seq, p_act, desc, score

popN = cfg.ctcm.population;
T0 = cfg.ctcm.T0;
Tmin = cfg.ctcm.Tmin;
H = cfg.ctcm.anneal_horizon;
burnin = floor(cfg.ctcm.burnin_frac * cfg.ctcm.proposals_per_run);

% Initialize population by sampling from seeds + random feasible repairs
pop = strings(popN,1);
for i=1:popN
    s = seed_seqs(randi(numel(seed_seqs)));
    s = random_perturb(s, cfg, 3); % small random changes
    s = repair_to_windows(s, cfg, 2);
    pop(i) = s;
end

% Track bests
best_self = pop;
best_score = -inf(popN,1);

% tribe structure
numTr = cfg.ctcm.num_tribes;
tribeSize = popN / numTr;
if abs(tribeSize - round(tribeSize)) > 1e-6
    error('population must be divisible by num_tribes');
end
tribeSize = round(tribeSize);

archive = struct('seq',{},'p_act',{},'desc',{},'score',{});

global_best_seq = pop(1);
global_best_score = -inf;

for t = 1:cfg.ctcm.proposals_per_run
    % anneal temperature
    T = max(Tmin, T0 * exp(-3.0 * t / H));

    % coarse-to-fine mixing coefficient
    if t <= burnin
        mix = 0.0;
    else
        mix = (t - burnin) / max(1, (cfg.ctcm.proposals_per_run - burnin));
        mix = min(max(mix,0),1);
    end

    for i=1:popN
        % Identify tribe
        tr = ceil(i / tribeSize);
        idx0 = (tr-1)*tribeSize + 1;
        idx1 = tr*tribeSize;
        tribe_idxs = idx0:idx1;

        % Chief: best in tribe so far
        [~, rel] = max(best_score(tribe_idxs));
        chief_idx = tribe_idxs(rel);
        chief_seq = best_self(chief_idx);

        % Anchors mixture
        r = rand();
        if r < 0.34
            anchor = best_self(i);
        elseif r < 0.67
            anchor = chief_seq;
        else
            anchor = global_best_seq;
        end

        s = pop(i);

        % Propose edit biased toward anchor + repair
        s_prop = propose_toward_anchor(s, anchor, cfg);
        s_prop = repair_to_windows(s_prop, cfg, 2);

        % Score (coarse-to-fine)
        [R_prop, p_act, desc] = objective_score(s_prop, seed_seqs, model, cfg, mix);
        [R_cur, ~, ~] = objective_score(s, seed_seqs, model, cfg, mix);

        accept = true;
        if cfg.ctcm.metropolis
            if R_prop < R_cur
                accept = (rand() < exp((R_prop - R_cur)/max(1e-6,T)));
            end
        else
            accept = (R_prop >= R_cur);
        end

        if accept
            pop(i) = s_prop;
            if R_prop > best_score(i)
                best_score(i) = R_prop;
                best_self(i) = s_prop;
            end
            if R_prop > global_best_score
                global_best_score = R_prop;
                global_best_seq = s_prop;
            end

            % Final gate for archive admission (use calibrated p_act + conformal + windows)
            p_act_cal = sigmoid(predict_logit(model.activity, build_features(s_prop, cfg)) / model.temperature);
            % No hemolysis in this no-external core run
            ok = conformal_accept(p_act_cal, model.conformal, 0.0, cfg, desc);

            if ok
                archive = try_archive_add(archive, s_prop, p_act_cal, desc, R_prop, cfg);
            end
        end
    end
end

% Sort archive by score desc
if ~isempty(archive)
    [~,ord] = sort([archive.score], 'descend');
    archive = archive(ord);
end

end

function s = random_perturb(s, cfg, k)
for i=1:k
    s = mutate_one(s, cfg);
end
end

function s = mutate_one(s, cfg)
% single random edit following p_sub/p_ins/p_del
r = rand();
if r < cfg.ctcm.p_sub
    pos = randi(strlength(s));
    aa = random_aa();
    s = replaceBetween(s, pos, pos, aa);
elseif r < cfg.ctcm.p_sub + cfg.ctcm.p_ins
    if strlength(s) >= cfg.data.max_len, return; end
    pos = randi(strlength(s)+1);
    aa = random_aa();
    s = insertBefore(s, pos, aa);
else
    if strlength(s) <= cfg.data.min_len, return; end
    pos = randi(strlength(s));
    s = eraseBetween(s, pos, pos);
end
end

function aa = random_aa()
aa20 = 'ACDEFGHIKLMNPQRSTVWY'; % char
aa = string(aa20(randi(numel(aa20))));
end
function s2 = propose_toward_anchor(s, anchor, cfg)
% Choose an edit type and bias substitutions toward matching anchor at mismatched positions.
s2 = s;

% If lengths differ, sometimes insert/delete to approach anchor length
if strlength(s) ~= strlength(anchor) && rand() < 0.4
    if strlength(s) < strlength(anchor) && strlength(s) < cfg.data.max_len
        pos = randi(strlength(s)+1);
        % insert residue sampled from anchor context if possible
        apos = min(pos, strlength(anchor));
        aa = extractBetween(anchor, apos, apos);
        if aa == ""; aa = random_aa(); end
        s2 = insertBefore(s2, pos, aa);
        return;
    elseif strlength(s) > strlength(anchor) && strlength(s) > cfg.data.min_len
        pos = randi(strlength(s));
        s2 = eraseBetween(s2, pos, pos);
        return;
    end
end

% Otherwise substitute at a mismatched position if possible
L = strlength(s);
L2 = strlength(anchor);
Lmin = min(L, L2);
mismatch = [];
for i=1:Lmin
    if extractBetween(s,i,i) ~= extractBetween(anchor,i,i)
        mismatch(end+1) = i; %#ok<AGROW>
    end
end
if ~isempty(mismatch)
    pos = mismatch(randi(numel(mismatch)));
    aa = extractBetween(anchor, pos, pos);
else
    pos = randi(L);
    aa = random_aa();
end
s2 = replaceBetween(s2, pos, pos, aa);

% Occasionally apply a second random mutation to keep exploration
if rand() < 0.2
    s2 = mutate_one(s2, cfg);
end
end

function s = repair_to_windows(s, cfg, steps)
% Greedy repair steps to move into descriptor windows.
for k=1:steps
    d = compute_physchem(s);
    if d.length < cfg.data.min_len
        s = insertBefore(s, 1, "K"); continue;
    end
    if d.length > cfg.data.max_len
        s = eraseBetween(s, strlength(s), strlength(s)); continue;
    end

    % Charge
    if d.charge < cfg.windows.charge_min
        s = swap_to_more_cationic(s); continue;
    end
    if d.charge > cfg.windows.charge_max
        s = swap_to_less_cationic(s); continue;
    end

    % GRAVY
    if d.gravy < cfg.windows.gravy_min
        s = swap_to_more_hydrophobic(s); continue;
    end
    if d.gravy > cfg.windows.gravy_max
        s = swap_to_less_hydrophobic(s); continue;
    end

    % Boman
    if d.boman > cfg.windows.boman_max
        s = swap_to_reduce_boman(s); continue;
    end
end
end

function s = swap_to_more_cationic(s)
pos = randi(strlength(s));
opts = 'KR'; aa = string(opts(randi(numel(opts))));
s = replaceBetween(s, pos, pos, aa);
end

function s = swap_to_less_cationic(s)
pos = randi(strlength(s));
opts = 'ADEGNPQSTV'; aa = string(opts(randi(numel(opts))));
s = replaceBetween(s, pos, pos, aa);
end

function s = swap_to_more_hydrophobic(s)
pos = randi(strlength(s));
opts = 'AILMFWV'; aa = string(opts(randi(numel(opts))));
s = replaceBetween(s, pos, pos, aa);
end

function s = swap_to_less_hydrophobic(s)
pos = randi(strlength(s));
opts = 'DEHKNPQRST'; aa = string(opts(randi(numel(opts))));
s = replaceBetween(s, pos, pos, aa);
end

function s = swap_to_reduce_boman(s)
% Heuristic: replace random residue with hydrophobic to lower solubility average
pos = randi(strlength(s));
opts = 'AILVFM'; aa = string(opts(randi(numel(opts))));
s = replaceBetween(s, pos, pos, aa);
end

function [R, p_act, desc] = objective_score(seq, ref_seqs, model, cfg, mix)
% mix=0: exploration reward (novelty - penalties)
% mix=1: add calibrated activity evidence

desc = compute_physchem(seq);

% penalty for leaving windows (smooth hinge-like quadratic)
pen = 0;
pen = pen + quad_hinge(desc.charge, cfg.windows.charge_min, cfg.windows.charge_max);
pen = pen + quad_hinge(desc.gravy,  cfg.windows.gravy_min,  cfg.windows.gravy_max);
pen = pen + quad_hinge(desc.boman,  -inf,                 cfg.windows.boman_max);
% length regularizer around mean prior (use 20)
len0 = 20;
pen = pen + cfg.ctcm.lambda_len * (desc.length - len0).^2;

% novelty: min distance in embedding space to reference set (seeds)
nov = novelty_distance(seq, ref_seqs, cfg);

% calibrated activity probability (optional)
Xq = build_features(string(seq), cfg);
z = predict_logit(model.activity, Xq);
p_act = sigmoid(z / model.temperature);

R_explore = cfg.ctcm.lambda_novel * nov - cfg.ctcm.lambda_penalty * pen;
R_refine  = R_explore + 2.0 * p_act; % activity bonus weight (tuneable)

R = (1-mix)*R_explore + mix*R_refine;
end

function v = quad_hinge(x, lo, hi)
if isfinite(lo) && x < lo
    v = (lo - x).^2;
elseif isfinite(hi) && x > hi
    v = (x - hi).^2;
else
    v = 0;
end
end

function d = novelty_distance(seq, ref_seqs, cfg)
% Use embedding distance (same as model embedding for reproducibility)
e = seq_embedding(string(seq), cfg);
if isempty(e)
    d = 0;
    return;
end
mind = inf;
for i=1:numel(ref_seqs)
    e2 = seq_embedding(ref_seqs(i), cfg);
    dd = norm(e - e2, 2);
    if dd < mind
        mind = dd;
    end
end
% Convert to a reward-like metric (bigger is more novel)
d = mind;
end

function archive = try_archive_add(archive, seq, p_act, desc, score, cfg)
% Add if archive not full or better than worst; prune near duplicates.
% Near-duplicate check: Hamming identity if same length.

% prune duplicates vs archive
for i=1:numel(archive)
    if strlength(seq) == strlength(archive(i).seq)
        id = mean(char(seq) == char(archive(i).seq));
        if id >= cfg.ctcm.min_hamming_identity
            return;
        end
    end
end

cand = struct('seq', string(seq), 'p_act', p_act, 'desc', desc, 'score', score);

if numel(archive) < cfg.ctcm.archive
    archive(end+1) = cand;
else
    [worstScore, wi] = min([archive.score]);
    if score > worstScore
        archive(wi) = cand;
    end
end
end
