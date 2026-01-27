function X = build_features(seqs, cfg)
%BUILD_FEATURES Build feature matrix from sequences.
% X = [physchem_zscored , embedding]

n = numel(seqs);
phys = zeros(n,4); % [len, charge, gravy, boman]
emb_cell = cell(n,1);
emb_dim = [];

for i=1:n
    d = compute_physchem(seqs(i));
    phys(i,:) = [d.length, d.charge, d.gravy, d.boman];
    e = seq_embedding(seqs(i), cfg);
    emb_cell{i} = e;
    if isempty(emb_dim); emb_dim = numel(e); end
end

if isempty(emb_dim); emb = zeros(n,0);
else
    emb = zeros(n, emb_dim);
    for i=1:n
        emb(i,:) = emb_cell{i};
    end
end

if cfg.model.zscore
    mu = mean(phys,1);
    sig = std(phys,0,1) + 1e-12;
    phys = (phys - mu) ./ sig;
end

X = [phys, emb];
end
