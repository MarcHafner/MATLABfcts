function mat_out = NormSumUnit(mat_in, dim)
% mat_out = NormSumUnit(mat_in, dim)
%   normalized each column (or along dimension  dim  ) to have a sum of 1
%

if ~exist('dim','var')
    dim = 1;
end
dims = ones(1,length(size(mat_in)));
dims(dim) = size(mat_in,dim);

mat_out = mat_in ./ repmat(sum(mat_in,dim), dims);
