function t_fits = extract_IC50(t_data, keys, pcutoff)

%%
if ~exist('keys','var')
    keys = intersect(t_data.Properties.VariableNames, ...
        {'CellLine' 'DrugName' 'Time' 'SeedingNumber' 'date'});
end

if ~exist('pcutoff','var')
    pcutoff = .05;
end
fitopt.pcutoff = pcutoff;


t_keys = unique(t_data(:,keys));

DoGI50 = ismember('RelGrowth', t_data.Properties.VariableNames);

t_fits = table;
for ik = 1:height(t_keys)
    loop_waitbar(ik, height(t_keys))
    %%
    subt = t_data(eqtable(t_keys(ik,:), t_data(:,keys)),:);


    [IC50, Hill, Emax, Area, r2, EC50, fit] = ...
        ICcurve_fit(subt.Conc, subt.RelCellCnt, 'IC50', fitopt);

    t_temp = [t_keys(ik,:) table(IC50, Hill, Emax, Area, r2, EC50) ...
        table({fit}, {subt.Conc'}, {subt.RelCellCnt'}, 'VariableNames', ...
        {'fit' 'Conc' 'RelCellCnt'})];

    if DoGI50
        [GI50, ~, GImax, GIArea, GI_r2, ~, GI_fit] = ...
            ICcurve_fit(subt.Conc, subt.RelGrowth, 'GI50', fitopt);

        t_temp = [t_temp table(GI50, GImax, GIArea, GI_r2) ...
        table({GI_fit}, {subt.RelGrowth'}, 'VariableNames', {'GI_fit' 'RelGrowth'})];

    end

    t_fits = [t_fits; t_temp];

end
