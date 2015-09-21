function [n, cats] = hist_cat(cat_array, cats)

if ~exist('cats','var')
    cats = unique(cat_array);
end
n = NaN(length(cats),size(cat_array,2));
for i = 1:length(n)
    n(i,:) = sum(cat_array==cats(i));
end
