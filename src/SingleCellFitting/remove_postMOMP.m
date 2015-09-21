function compiled_data = remove_postMOMP(compiled_data)

lD = height(compiled_data.ExpKey);
fields = setdiff(fieldnames(compiled_data.Traj), 'T')';

for iE = 1:lD

    filter = ~compiled_data.Fate(iE).Surviving;

    for f = fields
        if all(size(compiled_data.Traj(iE).(f{:})) == [length(compiled_data.Traj(iE).T), length(filter)])
            compiled_data.Traj(iE).([f{:} '_postMOMP'])=...
                compiled_data.Traj(iE).(f{:});
        end
        for j=find(filter)
            idx = compiled_data.Fate(iE).FRETPreMompTimeIdx(j);

            if isinf(idx)
                continue
            end
            compiled_data.Traj(iE).([f{:} '_postMOMP'])(1:idx,j)=NaN;
            compiled_data.Traj(iE).(f{:})((idx+1):end,j)=NaN;
        end
    end
end
