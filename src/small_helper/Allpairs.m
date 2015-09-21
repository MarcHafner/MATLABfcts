% pairs = Allpairs(vec1, vec2)
%   e.g. vec1=[1 2], vec2=[4 5];
%       -> pairs [1 4; 1 5; 2 4; 2 5]
%
function pairs = Allpairs(vec1, vec2)

if nargin==1
    idx = size(vec1)==2;
    assert(length(idx)==2 && any(idx), 'Need a 2 column matrix or 2 vectors')
    if sum(idx)==1
        if idx(1)
            vec2 = vec1(2,:);
            vec1 = vec1(1,:);
        else
            vec2 = vec1(:,1);
            vec1 = vec1(:,2);
        end
    end
else
    assert(isvector(vec1) && isvector(vec2), 'Need a 2 column matrix or 2 vectors')
end


[p,q] = meshgrid(vec1, vec2);
pairs = [p(:) q(:)];
