function CV = CoeffVar(x, dim)

if ~exist('dim','var')
    dim = 1;
end

CV = std(x,[],dim)./mean(x,dim);
