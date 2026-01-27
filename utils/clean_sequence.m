function s = clean_sequence(s, cfg)
% Clean whitespace and enforce canonical AA.
s = upper(string(s));
s = erase(s, " ");
s = erase(s, char(9)); % tabs
s = erase(s, char(10));
s = erase(s, char(13));

if cfg.data.allow_ambiguous
    return;
end
% Keep only canonical 20 AA
mask = regexp(s, "^[ACDEFGHIKLMNPQRSTVWY]+$");
if isempty(mask)
    s = "";
end
end
