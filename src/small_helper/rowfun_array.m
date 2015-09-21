function b = rowfun_array(fct, a)
% b = rowfun_array(fct, a)
%   apply the function fct on each row of the matrix a to yield a column
%   vector b

b = NaN(size(a,1),1);
for i=1:length(b)
    b(i) = fct(a(i,:));
end
