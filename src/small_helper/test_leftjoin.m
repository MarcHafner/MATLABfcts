function [t_out, Ridx, idx] = test_leftjoin(t_left, t_right, varargin)
% [t_out, Ridx, idx] = leftjoin(t_left, t_right, varargin)
%   outerjoin for talbe to the left, maintain order of t_left
%
%   Ridx are the indexes in t_right corresponding to the rows of t_left
%

[t_out, idx, Ridx] = outerjoin(t_left, t_right, 'type', 'left', ...
    'MergeKeys', true, varargin{:});

t_out = t_out(sortidx(idx),:);
Ridx = Ridx(sortidx(idx));

nullidx = histcounts(idx,1:height(t_left))==0;
if any(nullidx)
    disp('Rows in left table Missing right term:')
    t_left(nullidx,:)
end


multiidx = find(histcounts(idx,1:height(t_left))>1);
if any(multiidx)
    disp('Rows in left table with multiple right terms:')
    for i=1:length(multiidx)
        fprintf('---   %i    -----\n', i)
        disp(t_left(multiidx(i),:))
        disp(t_right(Ridx(multiidx(i)==idx),:))

        if mod(i,10)==0
            pause
        end
    end
end
