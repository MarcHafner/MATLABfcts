function tab = CheckVariableNames(tab, default_VariableNames, named_VariableNames)
% tab = CheckVariableNames(tab, default_VariableNames, named_VariableNames)
%
%   check a table for the variable names and correct the
%   'named_VariableNames'
%       tab is a table
%       default_VariableNames is a list of variable names
%       named_VariableNames (optional) is a n x 2 list of variable names (1st column
%           is the original names, 2nd is the changed names)
%
%

for i=1:length(default_VariableNames)
    assert(ismember(default_VariableNames{i}, tab.Properties.VariableNames), ...
        ['columns ' default_VariableNames{i} ' missing in table'])
end

if exist('named_VariableNames','var')
    for i=1:size(named_VariableNames,1)
        assert(ismember(named_VariableNames{i,1}, tab.Properties.VariableNames), ...
            ['columns ' named_VariableNames{i,1} ' missing in table'])
        tab.Properties.VariableNames{named_VariableNames{i,1}} = ...
            named_VariableNames{i,2};
    end
end
