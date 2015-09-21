function output = arrayfun_brdcast(fct, varargin)
% output = arrayfun_brdcast(fct, varargin)
%   rewritten arrayfun to broadcast any single cell input
%

assert(isa(fct, 'function_handle'))

for i=1:length(varargin)
    if isa(varargin{i},'table')
        varargin{i} = table2cell(varargin{i});
    end
end

varargin2 = varargin;
maxlength = max( cellfun(@length,varargin) .* (cellfun(@iscell,varargin)));
for i=1:length(varargin)
    if ~iscell(varargin{i}) || length(varargin{i})~=maxlength
        assert( ~iscell(varargin{i}) || length(varargin{i})==1 )

        if isvector(varargin{i}) && length(varargin{i})==maxlength
            varargin2{i} = mat2cell(varargin{i});
        else
            varargin2{i} = cell(maxlength,1);
            for j=1:maxlength
                varargin2{i}{j} = varargin{i};
            end
        end
    elseif isrow(varargin{i})
        varargin2{i} = varargin{i}';
    end
end

output = arrayfun(fct, varargin2{:}, 'uniformoutput', false);
