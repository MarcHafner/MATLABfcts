function b = rowfun_cells(fct, a)
% b = rowfun_cells(fct, a)
%   apply the function fct on each row of the cell array a to yield a column
%   cell vector b

b = cell(size(a,1),1);
for i=1:length(b)
    b{i} = fct(a(i,:));
end
