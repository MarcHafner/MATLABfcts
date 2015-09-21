function mat = struct2mat(s, field)

if isnumeric(s(1).(field))
    mat = NaN(size(s));
    for i=1:numel(s)
        if ~isempty(s(i).(field))
            mat(i) = s(i).(field);
        end
    end
else

    mat = cell(size(s));
    for i=1:numel(s)
        if ~isempty(s(i).(field))
            mat{i} = s(i).(field);
        end
    end
end
