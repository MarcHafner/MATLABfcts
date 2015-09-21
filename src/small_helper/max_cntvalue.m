function maxval = max_cntvalue(a, dim)
% maxval = max_cntvalue(a, dim)
%   fund the value with the most occurences in each column of a

if exist('dim','var') && dim==2
    a = a';
end

maxval = NaN(1, size(a,2));
for i = 1:size(a,2)
    vals = unique(a(:,i));
    if length(vals)==1
        maxval(i)=vals; continue
    end
    n = hist(a(:,i), vals);
    if sum(n==max(n))>1
        warnprintf('Multiple values with max occurence')
    end
    maxval(i) = vals(argmax(n));
end
