function name = MATLABsafename(name)
% name = MATLABsafename(name)
%   add x_ infront of a name starting with a number character
%

if ischar(name)
    if ~isempty(str2num(name(1)))
        name = ['x_' name];
    end
elseif iscellstr(name)
    for i=1:length(name)
        name{i} = MATLABsafename(name{i});
    end
end
