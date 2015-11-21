function  test_innerjoin(t_left, t_right, varargin)
% test_innerjoin(t_left, t_right, varargin)
%   look for missing matches or double matches
%
%

[t_out, Lidx, Ridx] = outerjoin(t_left, t_right, ...
    'MergeKeys', true, varargin{:});

nullidx = histcounts(Ridx,1:height(t_out))==0;
if any(nullidx)
    disp('Rows in left table missing right equivalence:')
    t_out(nullidx,:)
end


multiidx = find(histcounts(Lidx,1:height(t_out))>1);
if any(multiidx)
    disp('Rows in left table with multiple right terms:')
    t_out(multiidx,:)
end

nullidx = histcounts(Lidx,1:height(t_out))==0;
if any(nullLidx)
    disp('Rows in right table missing left equivalence:')
    t_out(nullidx,:)
end


multiidx = find(histcounts(Ridx,1:height(t_out))>1);
if any(multiidx)
    disp('Rows in right table with multiple left terms:')
    t_out(multiidx,:)
end
