function m = infmean(x,dim)
%INFMEAN Mean value, ignoring Infs.
%   M = INFMEAN(X) returns the sample mean of X, treating Infs as missing
%   values.  For vector input, M is the mean value of the non-Inf elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-Inf elements in each column.  For N-D arrays, INFMEAN operates
%   along the first non-singleton dimension.
%
%   INFMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, INFMEDIAN, INFSTD, INFVAR, INFMIN, INFMAX, INFSUM.

%   Copyright 1993-2004 The MathWorks, Inc.


% Find Infs and set them to zero
if nargin == 1
    minf = mean(x);
else
    minf = mean(x,dim);
end

nans = isinf(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-Infs.
    n = sum(~nans);
    n(n==0) = Inf; % prevent divideByZero warnings
    % Sum up non-Infs, and divide by the number of non-Infs.
    m = sum(x) ./ n;
else
    % Count up non-Infs.
    n = sum(~nans,dim);
    n(n==0) = Inf; % prevent divideByZero warnings
    % Sum up non-Infs, and divide by the number of non-Infs.
    m = sum(x,dim) ./ n;
end

m(n==0) = minf(n==0); % conserve vectors with all Inf or all -Inf
