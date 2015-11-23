function mat_out = NormNormUnit(mat_in, dim)
% mat_out = NormNormUnit(mat_in, dim)
%   normalized each column (or along dimension  dim  ) to have a L2-norm of 1
%

if ~exist('dim','var')
    dim = 1;
end
dims = ones(1,length(size(mat_in)));
dims(dim) = size(mat_in,dim);

mat_out = mat_in ./ repmat(sqrt(sum(mat_in.^2,dim)), dims);
