function str = AnyToString(in)

if iscell(in)
    str = cellfun2(@AnyToString, in);
    return
elseif iscategorical(in)
    str = cell(size(in));
    for i=1:numel(in)
        str{i} = char(in(i));
    end
elseif isnumeric(in)
    str = cellfun2(@num2str, num2cell(in));
elseif ischar(in)
    str = in;
    return
else
    str = evalc('disp(in)');
end

if iscell(str) && isscalar(str)
    str = str{:};
end
