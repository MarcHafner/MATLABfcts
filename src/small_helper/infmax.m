function [varargout] = infmax(x, varargin)
%INFMAX Maximum value, ignoring NaNs.
%   M = INFMAX(A) returns the maximum of A with NaNs treated as missing.
%   For vectors, M is the largest non-NaN element in A.  For matrices, M is
%   a row vector containing the maximum non-NaN element from each column.
%   For N-D arrays, INFMAX operates along the first non-singleton
%   dimension.
%
%   [M,NDX] = INFMAX(A) returns the indices of the maximum values in A.  If
%   the values along the first non-singleton dimension contain more than
%   one maximal element, the index of the first one is returned.
%
%   M = INFMAX(A,B) returns an array the same size as A and B with the
%   largest elements taken from A or B.  Either one can be a scalar.
%
%   [M,NDX] = INFMAX(A,[],DIM) operates along the dimension DIM.
%
%   See also MAX, INFMIN, INFMEAN, INFMEDIAN, INFMIN, INFVAR, INFSTD.

%   Copyright 1993-2004 The MathWorks, Inc.


% Call [m,ndx] = max(a,b) with as many inputs and outputs as needed
x(isinf(x)) = nan;
[varargout{1:nargout}] = max(x, varargin{:});
