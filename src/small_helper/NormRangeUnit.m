function mat_out = NormRangeUnit(mat_in, dim)
% mat_out = NormRangeUnit(mat_in, dim)
%   normalized each column (or along dimension  dim  ) to range from 0 to 1
%

if ~exist('dim','var')
    dim = 1;
end
dims = ones(1,length(size(mat_in)));
dims(dim) = size(mat_in,dim);

mat_in = mat_in - repmat(min(mat_in,[],dim), dims);
mat_out = mat_in ./ repmat(max(mat_in,[],dim), dims);
