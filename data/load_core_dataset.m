function [seqs, y, meta] = load_core_dataset(pos_xlsx, neg_xlsx, cfg)
%LOAD_CORE_DATASET Load positives and unlabeled/background as core corpus.
% pos_xlsx: expects a column containing sequences (preferably named 'Sequence')
% neg_xlsx: expects a column containing sequences (preferably named 'Sequence')
%
% Important: the paper constrains lengths to [min_len,max_len]. In your positive file
% some sequences are > max_len. To keep the provided positive count while enforcing
% length<=max_len, cfg.data.long_strategy controls what happens.

Tpos = readtable(pos_xlsx);
Tneg = readtable(neg_xlsx);

seq_pos = pick_sequence_column(Tpos);
seq_neg = pick_sequence_column(Tneg);

% Clean (canonical AA filter)
seq_pos_raw = string(seq_pos);
seq_neg_raw = string(seq_neg);

seq_pos_clean = arrayfun(@(x) clean_sequence(x, cfg), seq_pos_raw);
seq_neg_clean = arrayfun(@(x) clean_sequence(x, cfg), seq_neg_raw);

% Handle long sequences (positives + negatives) before length filtering
seq_pos_clean = handle_long_sequences(seq_pos_clean, cfg);
seq_neg_clean = handle_long_sequences(seq_neg_clean, cfg);

% Filter by length and drop empties
pos_ok = (seq_pos_clean ~= "") & (strlength(seq_pos_clean) >= cfg.data.min_len) & (strlength(seq_pos_clean) <= cfg.data.max_len);
neg_ok = (seq_neg_clean ~= "") & (strlength(seq_neg_clean) >= cfg.data.min_len) & (strlength(seq_neg_clean) <= cfg.data.max_len);

seq_pos = seq_pos_clean(pos_ok);
seq_neg = seq_neg_clean(neg_ok);

seqs = [seq_pos; seq_neg];
y = [ones(numel(seq_pos),1); zeros(numel(seq_neg),1)];

% Deduplication (exact duplicates)
before_dedup = numel(y);
if cfg.data.dedup
    [seqs, ia] = unique(seqs, 'stable');
    y = y(ia);
end
after_dedup = numel(y);

meta = struct();
meta.n_pos = sum(y==1);
meta.n_u   = sum(y==0);
meta.total = numel(y);

% Report attrition
meta.input_pos = height(Tpos);
meta.input_u = height(Tneg);
meta.pos_noncanonical_or_empty = sum(seq_pos_clean=="");
meta.u_noncanonical_or_empty   = sum(seq_neg_clean=="");
meta.pos_dropped_len = sum(~pos_ok & seq_pos_clean~="");
meta.u_dropped_len   = sum(~neg_ok & seq_neg_clean~="");
meta.dedup_removed   = before_dedup - after_dedup;
meta.long_strategy   = cfg.data.long_strategy;

fprintf('Loaded core dataset: P=%d, U=%d, total=%d (long_strategy=%s)\n', meta.n_pos, meta.n_u, meta.total, cfg.data.long_strategy);

% Persist a report (useful when positives are unexpectedly reduced)
if ~isfolder(cfg.out_dir); mkdir(cfg.out_dir); end
fid=fopen(fullfile(cfg.out_dir,'data_report.json'),'w');
fprintf(fid, '%s', jsonencode(meta));
fclose(fid);

end

function col = pick_sequence_column(T)
vars = string(T.Properties.VariableNames);
idx = find(vars=="Sequence", 1);
if isempty(idx)
    % fallback: first column
    idx = 1;
end
col = T{:,idx};
end

function seqs = handle_long_sequences(seqs, cfg)
seqs = string(seqs);
maxL = cfg.data.max_len;
if cfg.data.long_strategy == "drop"
    % do nothing; length filter later will remove
    return;
elseif cfg.data.long_strategy == "truncate"
    for i=1:numel(seqs)
        if strlength(seqs(i)) > maxL
            seqs(i) = extractBetween(seqs(i), 1, maxL);
        end
    end
elseif cfg.data.long_strategy == "best_window"
    for i=1:numel(seqs)
        if strlength(seqs(i)) > maxL
            seqs(i) = best_window_crop(seqs(i), maxL);
        end
    end
else
    error('Unknown cfg.data.long_strategy=%s', cfg.data.long_strategy);
end
end

function w = best_window_crop(s, Lw)
% Choose the Lw-length window maximizing a simple AMP-like heuristic
% score = charge - 0.5*abs(gravy) - 0.5*boman
s = string(s);
L = strlength(s);
bestScore = -inf;
best = extractBetween(s, 1, Lw);
for start = 1:(L - Lw + 1)
    cand = extractBetween(s, start, start+Lw-1);
    d = compute_physchem(cand);
    score = d.charge - 0.5*abs(d.gravy) - 0.5*d.boman;
    if score > bestScore
        bestScore = score;
        best = cand;
    end
end
w = best;
end
