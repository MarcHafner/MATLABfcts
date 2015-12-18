function out = peek(x, n)
% out = peek(x, n)
%   relies on the function 'head' from datarail

if nargin==1
    n = 10;
end

[l, d] = length_(x);
idx = unique([1:max(1,floor(l/n)):l l]);

if istable(x)
    out = x(idx, :);
else
    sz = size(x);
    sz(d) = n;
    out = reshape(x, l, []);
    out = reshape(out(idx, :), sz);
end
