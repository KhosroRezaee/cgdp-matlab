function splits = make_cluster_split(seqs, y, cfg)
%MAKE_CLUSTER_SPLIT Cluster split at <=40% identity.
% Paper-faithful mode uses CD-HIT. If CD-HIT is unavailable, fall back to a
% MATLAB-only greedy clustering (approximate; may differ from CD-HIT).

rng(cfg.cluster.seed);

if ~isfolder(cfg.out_dir); mkdir(cfg.out_dir); end
fasta = export_fasta(seqs, fullfile(cfg.out_dir, "core.fasta"));

clusters = {};

use_cdhit = cfg.cluster.use_cdhit;
if use_cdhit
    % Quick availability check
    if ~command_exists(cfg.cluster.cdhit_bin)
        warning('CD-HIT not found on PATH. Falling back to MATLAB clustering.');
        use_cdhit = false;
    end
end

if use_cdhit
    out_prefix = fullfile(cfg.out_dir, "cdhit40");
    cmd = sprintf('%s -i "%s" -o "%s" -c %.2f -n 2 -d 0 -T 0 -M 0', ...
        cfg.cluster.cdhit_bin, fasta, out_prefix, cfg.cluster.cdhit_identity);
    st = system(cmd);
    if st ~= 0
        warning('CD-HIT failed; falling back to MATLAB clustering. Command was: %s', cmd);
        use_cdhit = false;
    else
        clstr = out_prefix + ".clstr";
        clusters = parse_cdhit_clstr(clstr);
    end
end

if ~use_cdhit
    clusters = greedy_cluster_identity(seqs, cfg.cluster.cdhit_identity);
end

% Shuffle clusters then allocate by fraction
m = numel(clusters);
ord = randperm(m);

n_tr = max(1, round(cfg.cluster.train_frac * m));
n_va = max(1, round(cfg.cluster.val_frac * m));
n_te = m - n_tr - n_va;
if n_te < 1
    n_te = 1; n_va = max(1, m-n_tr-1);
end

tr_clusters = ord(1:n_tr);
va_clusters = ord(n_tr+1:n_tr+n_va);
te_clusters = ord(n_tr+n_va+1:end);

train_idx = flatten_clusters(clusters, tr_clusters);
val_idx   = flatten_clusters(clusters, va_clusters);
test_idx  = flatten_clusters(clusters, te_clusters);

% Enforce minimum positives in val/test by moving whole clusters from train if needed
[tr_clusters, va_clusters, te_clusters] = enforce_min_pos(clusters, y, tr_clusters, va_clusters, te_clusters, cfg);
train_idx = flatten_clusters(clusters, tr_clusters);
val_idx   = flatten_clusters(clusters, va_clusters);
test_idx  = flatten_clusters(clusters, te_clusters);
splits = struct();
splits.clusters = clusters;
splits.train_idx = sort(train_idx(:));
splits.val_idx   = sort(val_idx(:));
splits.test_idx  = sort(test_idx(:));

splits.n_train = numel(splits.train_idx);
splits.n_val   = numel(splits.val_idx);
splits.n_test  = numel(splits.test_idx);

splits.pos_train = sum(y(splits.train_idx)==1);
splits.pos_val   = sum(y(splits.val_idx)==1);
splits.pos_test  = sum(y(splits.test_idx)==1);

splits.used_cdhit = use_cdhit;

fprintf('Split sizes: train=%d (P=%d)  val=%d (P=%d)  test=%d (P=%d)  use_cdhit=%d\n', ...
    splits.n_train, splits.pos_train, splits.n_val, splits.pos_val, splits.n_test, splits.pos_test, splits.used_cdhit);
end

function clusters = greedy_cluster_identity(seqs, thr)
% Greedy clustering approximating CD-HIT at identity threshold thr.
% Representative = first sequence of cluster.
%
% Identity computed as LCS_length / max(len1,len2) (conservative).

seqs = string(seqs);
n = numel(seqs);

% Sort by length desc (CD-HIT-style)
lens = strlength(seqs);
[~, ord] = sort(lens, 'descend');

reps = strings(0,1);
clusters = {};

for kk = 1:n
    i = ord(kk);
    s = seqs(i);

    assigned = false;
    for c = 1:numel(reps)
        id = lcs_identity(char(s), char(reps(c)));
        if id >= thr
            clusters{c} = [clusters{c}(:); i]; %#ok<AGROW>
            assigned = true;
            break;
        end
    end

    if ~assigned
        reps(end+1) = s; %#ok<AGROW>
        clusters{end+1} = i; %#ok<AGROW>
        clusters{end} = clusters{end}(:);
    end
end
end

function id = lcs_identity(a, b)
% LCS-based identity proxy
la = length(a); lb = length(b);
L = lcs_len(a,b);
id = L / max(la, lb);
end

function L = lcs_len(a,b)
la=length(a); lb=length(b);
dp = zeros(lb+1,1);
for i=1:la
    prev = 0;
    for j=1:lb
        tmp = dp(j+1);
        if a(i)==b(j)
            dp(j+1) = prev + 1;
        else
            dp(j+1) = max(dp(j+1), dp(j));
        end
        prev = tmp;
    end
end
L = dp(lb+1);
end


function idx = flatten_clusters(clusters, ids)
% Flatten variable-length numeric vectors stored in clusters{...} into a single column.
parts = clusters(ids);
parts = cellfun(@(v) v(:), parts, 'UniformOutput', false);
idx = vertcat(parts{:});
end


function [tr, va, te] = enforce_min_pos(clusters, y, tr, va, te, cfg)
% Move clusters (whole) from train to val/test until minimum positives reached.
minV = cfg.cluster.min_pos_val;
minT = cfg.cluster.min_pos_test;

pos_in = @(ids) sum(y(flatten_clusters(clusters, ids))==1);

% Candidate clusters in train with positives, sorted by size ascending (small moves)
train_pos_clusters = tr(arrayfun(@(cid) sum(y(clusters{cid})==1)>0, tr));
sizes = arrayfun(@(cid) numel(clusters{cid}), train_pos_clusters);
[~,ord] = sort(sizes, 'ascend');
train_pos_clusters = train_pos_clusters(ord);

% Fill val
while pos_in(va) < minV && ~isempty(train_pos_clusters)
    cid = train_pos_clusters(1);
    train_pos_clusters(1) = [];
    tr(tr==cid) = [];
    va(end+1) = cid; %#ok<AGROW>
end

% Refresh remaining train positive clusters list
train_pos_clusters = tr(arrayfun(@(cid) sum(y(clusters{cid})==1)>0, tr));
sizes = arrayfun(@(cid) numel(clusters{cid}), train_pos_clusters);
[~,ord] = sort(sizes, 'ascend');
train_pos_clusters = train_pos_clusters(ord);

% Fill test
while pos_in(te) < minT && ~isempty(train_pos_clusters)
    cid = train_pos_clusters(1);
    train_pos_clusters(1) = [];
    tr(tr==cid) = [];
    te(end+1) = cid; %#ok<AGROW>
end
end
