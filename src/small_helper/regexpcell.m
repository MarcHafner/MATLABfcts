function pos = regexpcell(strcell, expression)
% pos = regexpcell(strcell,expression)
%   find the position of a expression expression in a cell of strings
%   strcell.
%   return the vector pos of the same length as strcell
%   containing the pos(i) of expression in strcell{i} or 0 if the
%   expression doesn't match strcell
%
%

pos = zeros(size(strcell));
temp = regexp(strcell, expression);
for i=1:length(temp)
    if ~isempty(temp{i})
        pos(i) = temp{i}(1);
    else
        pos(i) = 0;
    end
end
