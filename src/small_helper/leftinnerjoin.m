function [t_out, Lidx, Ridx] = leftinnerjoin(t_left, t_right, varargin)
% [t_out, Lidx, Ridx] = leftinnerjoin(t_left, t_right, varargin)
%
%

[t_out, Lidx, Ridx] = innerjoin(t_left, t_right, varargin{:});

t_out = t_out(sortidx(Lidx),:);
Ridx = Ridx(sortidx(Lidx));
Lidx = Lidx(sortidx(Lidx),:);

warnassert(all(hist(Lidx,1:height(t_left))<=1), ...
    'LEFTINNERJOIN: multiple entries from t_right match a row of t_left')

warnassert(all(hist(Ridx,1:height(t_right))<=1), ...
    'LEFTINNERJOIN: multiple entries from t_left match a row of t_right')
