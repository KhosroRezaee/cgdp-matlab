function fasta_path = export_fasta(seqs, out_path)
%EXPORT_FASTA Write sequences to FASTA file.
fid = fopen(out_path, 'w');
for i=1:numel(seqs)
    fprintf(fid, '>seq_%d\n%s\n', i, char(seqs(i)));
end
fclose(fid);
fasta_path = out_path;
end
