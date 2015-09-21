% [MI, pval] = MI_eval(x, y, nRand)
%
%   Calculating mutual information between two vectors x and y or a matrix
%   given in x; nRand is the number of permutation for evaluating the
%   p-value (default is 2000)
%


function [MI, pval] = MI_eval(x, y, nRand)

default_nRand = 2e3;
evalp = nargout>1;

% case of a single input that should be a matrix
if ~exist('y', 'var') || isempty(y) || (numel(y)==1 && max(size(x))>1)
    if exist('y', 'var') && numel(y)==1
        assert(exist('nRand','var'), ...
            'With x as a matrix, only one additional argument (nRand) can be given')
        nRand = y;
    elseif ~exist('nRand','var')
        nRand = default_nRand;
    end

    nDim = size(x,2);
    assert(nDim>1, 'If only x is given, it should have a least 2 columns');

    % this is valid as long as the bins contains the same number of points
    % (and not if they are equidistant)

    MI_1 = MI_eval(x(:,1), x(:,1));
    MI = zeros(nDim);

    if evalp;pval = zeros(nDim);end

    % test all pairwise combinations
    parfor i=1:nDim
        xi = x(:,i);
        t_MI = zeros(1,nDim); t_MI(i)=MI_1*.5;  % *0.5 is to correct for summing the transpose;
        t_pval = t_MI
        for j=(i+1):nDim
            if evalp    % evaluate the p-value
                [t_MI(j), t_pval(j)] = MI_eval(xi, x(:,j), nRand);
            else    % only the MI
                t_MI(j) = MI_eval(xi, x(:,j));
            end
        end
        if evalp;pval(i,:) = t_pval;end
        MI(i,:) = t_MI;

    end

    MI = MI+MI';
    if evalp; pval = pval+pval'; end

    return

else
    assert(length(x)==length(y), 'x and y should have the same length');
    x = ToColumn(x);
    y = ToColumn(y);
    if ~exist('nRand','var')
        nRand = default_nRand;
    end
end

total_pairs = length(x);

min_x = min(x);
min_y = min(y);

% number of bins used for discretization
% since this is a continuous measurement, need to bin to make it
% seem discreet
% metric_disc_number = 2*floor(log2(total_pairs)+1);
% metric_disc_number = ceil(sqrt(total_pairs));
metric_disc_number = max(4, ...
    min(floor(sqrt(total_pairs))-2, 2*floor(log2(total_pairs)+1)-2));

% finding adaptive bin edges            !!! if equally spaced, revise assumption for MI(x,x) !!!
% dist is the cumulative distribution value corresponding to x1_values
[dist1, x_vals] = ecdf(x);
[dist2, y_vals] = ecdf(y);

% creating bin edges with equal distribution of points in between
% must include minimum value of the array to be the left edge of first bin
x_bins = NaN(1, length(metric_disc_number));
y_bins = NaN(1, length(metric_disc_number));

for i=1:metric_disc_number
    x_bins(i) = x_vals(find(dist1 >= 1/metric_disc_number * i, 1, 'first'));
    y_bins(i) = y_vals(find(dist2 >= 1/metric_disc_number * i, 1, 'first'));
end

% such that the first value is slightly above the edge of the first bin
x_bins = [min_x - 1e-3; x_bins'];
y_bins = [min_y - 1e-3; y_bins'];

% not the cleanest to pre-assign the variables outside of the function but it
% is ~15% faster when evaluating the p-value
x_hist = NaN(length(x_bins)-1, length(x));
y_hist = x_hist;

MI = core_MI(x, y, x_bins, y_bins);

% getting p value for mutual information
% here a single p value is returned
if evalp

    % iterative version (faster)
    % --> not efficient to use a parfor here, too much overhead
    MI_postshuffle = NaN(1,nRand);
    for i=1:nRand
        MI_postshuffle(i) = core_MI(x, y(randperm(length(y))), x_bins, y_bins);
    end


    %% Real MI and Null MI comparison
    % # of times null hypothesis wins
    pval = sum(MI_postshuffle > MI) / nRand;

end



    function MI = core_MI(x, y, x_bins, y_bins)


        % x and y have the same number of bins
        for iBin=1:length(x_bins)-1
            x_hist(iBin, :) = (x > x_bins(iBin)) .* (x <= x_bins(iBin+1));
            y_hist(iBin, :) = (y > y_bins(iBin)) .* (y <= y_bins(iBin+1));
        end

        % executing matrix multiplication; order matters
        joint_prob = x_hist * y_hist' / length(x);

        % calculating marginal probabilities
        x_mprob = sum(joint_prob, 2);
        y_mprob = sum(joint_prob, 1);

        % calculating mutual information
        xy_mprobs_inv = x_mprob*y_mprob;

        idx = joint_prob>0;
        MI = sum(joint_prob(idx).*log2(joint_prob(idx)./xy_mprobs_inv(idx)));

    end

end

function col = ToColumn(vect)
% col = ToColumn(vect)
%   Convert a vector to a column vector

assert(any(size(vect)==1))
if isrow(vect)
    col=vect';
else
    col=vect;
end

end
