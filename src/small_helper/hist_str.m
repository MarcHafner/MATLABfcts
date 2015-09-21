function [n, unique_str] = hist_str(allstr, unique_str)
% [n, unique_str] = hist_str(allstr, unique_str)


if istable(allstr)
    assert(size(allstr,2)==1,'more than one variable in the table')
    allstr = allstr.(allstr.Properties.VariableNames{1});
end

if ~exist('unique_str', 'var')
    str_input = false;
    unique_str = unique(allstr);
else
    str_input = true;
end


[~,idx] = ismember(allstr, unique_str);

n = hist(idx, 0:length(unique_str));
n = n(2:end);
if ~str_input
    [n,order] = sort(n,'descend');
    unique_str = unique_str(order);
end
