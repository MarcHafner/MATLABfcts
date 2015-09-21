function [coeff, scores, explained, latent] = RankPCA(data, cutoff)
% [coeff, scores, explained, latent] = RankPCA(data, cutoff)
%   PCA analysis based on the gene ranks.
%   cutoff is used to filter the data with low standard deviation (default
%   value is 1e-6)
%


if ~exist('cutoff','var')
    cutoff = 1e-6;
end


coeff = zeros(size(data,2));
latent = zeros(size(data,2),1);

idx = std(data)>cutoff;

subdata = single(data(:,idx));
ranked = single(tiedrank(subdata));
clear subdata

[coeff(idx,1:sum(idx)),~, latent(1:sum(idx))] = princomp(ranked);
clear ranked

explained = 100*latent/sum(latent);

scores = (data-repmat(mean(data),size(data,1),1))*coeff;
