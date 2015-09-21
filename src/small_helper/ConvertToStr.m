function str = ConvertToStr(anything)

if iscell(anything) && numel(anything)>1
    str = cell(size(anything));
    for i=1:numel(anything)
        str{i} = ConvertToStr(anything{i});
    end
    return
elseif iscell(anything)
    str = ConvertToStr(anything{1});
    return
end

if isnumeric(anything)
    str = num2str(anything);
elseif iscategorical(anything)
    str=char(anything);
else
    str=anything;
end
