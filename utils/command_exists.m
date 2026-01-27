function tf = command_exists(cmd)
%COMMAND_EXISTS Check if a command exists on the system PATH.
% Works on Windows and Unix.
if ispc
    [st,~] = system(sprintf('where %s >NUL 2>NUL', cmd));
else
    [st,~] = system(sprintf('which %s >/dev/null 2>/dev/null', cmd));
end
tf = (st==0);
end
