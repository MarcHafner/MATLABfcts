function parts = regexpcellsplit(strcell, expression, selected)
% pos = regexpcellsplit(strcell, expression, selected)
%   find the parts in a cell of strings strcell split based on the regular
%   expression express. Return the cell array parts of the same length as
%   strcell containing the parts of expression in strcell{i} or the whole string
%   if the expression doesn't match strcell. If multiple
%   parts are found, all are passed unless selected is set to a integer
%   (e.g. 1 of only the first part of the split is reported). If selected
%   is a vector of integer, the selected parts are concatenated based on
%   the epxression
%
%   note: if selected is negative, works as -selected, BUT
%       and if -selected > length( splitted str), parts will be empty
%

parts = cell(size(strcell));
temp = regexp(strcell,expression,'split');
for i=1:numel(temp)
    if exist('selected','var')
        if length(selected)>1
            parts{i} = strjoin(temp{i}(unique(min(abs(selected), length(temp{i})))), ...
                expression);
        else
            if selected<0 && -selected>length(temp{i})
                parts{i} = '';
            else
                parts{i} = temp{i}{min(abs(selected), length(temp{i}))};
            end
        end
    elseif length(temp{i})==1
        parts{i} = temp{i}{1};
    else
        parts{i} = temp{i};
    end
end
