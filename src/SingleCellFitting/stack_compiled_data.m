function all_data = stack_compiled_data(subset_data1, subset_data2)
% all_data = stack_compiled_data(subset_data1, subset_data2)

assert(all(strcmp(fieldnames(subset_data1), fieldnames(subset_data2))))

all_data = subset_data1;
nExp = height(subset_data2.ExpKey);


var2to1 = setdiff(subset_data2.ExpKey.Properties.VariableNames, ...
    subset_data1.ExpKey.Properties.VariableNames);
if ~isempty(var2to1)
    all_data.ExpKey = [all_data.ExpKey array2table( NaN(height(all_data.ExpKey), ...
        length(var2to1)), 'VariableNames', var2to1)];
end

var1to2 = setdiff(subset_data1.ExpKey.Properties.VariableNames, ...
    subset_data2.ExpKey.Properties.VariableNames);
if ~isempty(var1to2)
    subset_data2.ExpKey = [subset_data2.ExpKey array2table( NaN(height(all_data.ExpKey), ...
        length(var1to2)), 'VariableNames', var1to2)];
end
subset_data2.ExpKey = subset_data2.ExpKey(:, all_data.ExpKey.Properties.VariableNames);

for f = fieldnames(subset_data2)'
    idx = find(size(subset_data2.(f{:}))==nExp,1,'first');
    if idx==1
        all_data.(f{:}) = [all_data.(f{:});
            subset_data2.(f{:})];
    elseif idx==2
        all_data.(f{:}) = [all_data.(f{:}) subset_data2.(f{:})];
    else
        error('Mismatched dimensions for field %s', f{:});
    end
end
