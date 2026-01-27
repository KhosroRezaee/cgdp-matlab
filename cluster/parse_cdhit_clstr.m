function clusters = parse_cdhit_clstr(clstr_path)
%PARSE_CDHIT_CLSTR Parse CD-HIT .clstr file into cell array of indices.
txt = fileread(clstr_path);
lines = splitlines(string(txt));
clusters = {};
cur = [];
for i=1:numel(lines)
    ln = strtrim(lines(i));
    if ln == ""
        continue;
    end
    if startsWith(ln, ">Cluster")
        if ~isempty(cur)
            clusters{end+1} = cur; %#ok<AGROW>
        end
        cur = [];
    else
        % Example: "0  30aa, >seq_12... *"
        tok = regexp(ln, ">seq_(\d+)", "tokens", "once");
        if ~isempty(tok)
            cur(end+1) = str2double(tok{1}); %#ok<AGROW>
        end
    end
end
if ~isempty(cur)
    clusters{end+1} = cur;
end
end
