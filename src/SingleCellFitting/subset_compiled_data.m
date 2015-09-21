function subset_data = subset_compiled_data(all_data, subsetIdx)
% subset_data = subset_compiled_data(all_data, subsetIdx)

nExp = height(all_data.ExpKey);

for f = fieldnames(all_data)'
    idx = find(size(all_data.(f{:}))==nExp,1,'first');
    if idx==1
        subset_data.(f{:}) = all_data.(f{:})(subsetIdx,:);
    elseif idx==2
        subset_data.(f{:}) = all_data.(f{:})(:,subsetIdx);
    else
        error('Mismatched dimensions for field %s', f{:});
    end
end
