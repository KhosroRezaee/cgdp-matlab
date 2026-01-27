function write_archive_csv(path, archive)
%WRITE_ARCHIVE_CSV Save archive to CSV
n = numel(archive);
seq = strings(n,1);
p  = zeros(n,1);
len = zeros(n,1);
chg = zeros(n,1);
gr  = zeros(n,1);
bo  = zeros(n,1);
sc  = zeros(n,1);

for i=1:n
    seq(i) = archive(i).seq;
    p(i) = archive(i).p_act;
    len(i) = archive(i).desc.length;
    chg(i) = archive(i).desc.charge;
    gr(i)  = archive(i).desc.gravy;
    bo(i)  = archive(i).desc.boman;
    sc(i)  = archive(i).score;
end

T = table(seq, p, sc, len, chg, gr, bo, ...
    'VariableNames', {'Sequence','P_active','Score','Length','Charge_pH7','GRAVY','BomanIndex'});
writetable(T, path);
end
