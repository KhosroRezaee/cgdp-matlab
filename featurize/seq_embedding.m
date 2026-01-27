function emb = seq_embedding(seq, cfg)
%SEQ_EMBEDDING Returns an embedding for a sequence.
% Paper uses frozen ProtT5 pooled states. Here we provide a local default:
%  - "dipeptide": AA composition (20) + dipeptide composition (400) = 420 dims
% Replace this function to load precomputed ProtT5 embeddings if desired.

seq = char(seq);
aa = 'ACDEFGHIKLMNPQRSTVWY';
L = length(seq);

switch cfg.model.embedding
    case "none"
        emb = [];
    case "dipeptide"
        % AA comp
        comp = zeros(1,20);
        for i=1:20
            comp(i) = sum(seq==aa(i)) / L;
        end
        % dipeptide comp
        di = zeros(1,400);
        if L >= 2
            for i=1:(L-1)
                a1 = find(aa==seq(i),1);
                a2 = find(aa==seq(i+1),1);
                if ~isempty(a1) && ~isempty(a2)
                    idx = (a1-1)*20 + a2;
                    di(idx) = di(idx) + 1;
                end
            end
            di = di / (L-1);
        end
        emb = [comp, di];
    otherwise
        error('Unknown cfg.model.embedding = %s', cfg.model.embedding);
end
end
