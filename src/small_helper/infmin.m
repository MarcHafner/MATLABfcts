function [varargout] = infmin(x, varargin)
%INFMIN Maximum value, ignoring NaNs.
%   M = INFMIN(A) returns the minimum of A with NaNs treated as missing.
%   For vectors, M is the largest non-NaN element in A.  For matrices, M is
%   a row vector containing the minimum non-NaN element from each column.
%   For N-D arrays, INFMIN operates along the first non-singleton
%   dimension.
%
%   [M,NDX] = INFMIN(A) returns the indices of the minimum values in A.  If
%   the values along the first non-singleton dimension contain more than
%   one minimal element, the index of the first one is returned.
%
%   M = INFMIN(A,B) returns an array the same size as A and B with the
%   largest elements taken from A or B.  Either one can be a scalar.
%
%   [M,NDX] = INFMIN(A,[],DIM) operates along the dimension DIM.
%
%   See also MIN, INFMIN, INFMEAN, INFMEDIAN, INFMIN, INFVAR, INFSTD.

%   Copyright 1993-2004 The MathWorks, Inc.


% Call [m,ndx] = min(a,b) with as many inputs and outputs as needed
x(isinf(x)) = nan;
[varargout{1:nargout}] = min(x, varargin{:});
