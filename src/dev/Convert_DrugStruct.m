function Convert_DrugStruct(file)

[pathstr,name,ext] = fileparts(file);

load(file)

%%

for i=1:length(drugs_struct)
    [DrugName, HMSLid] = splitHMSLid({drugs_struct(i).name});
    Drugs(i) = struct( ...
                'DrugName', DrugName, ...
                'HMSLid', HMSLid, ...
                'stock_conc', drugs_struct(i).nominal_conc, ...
                'layout', drugs_struct(i).layout );
end
Designs = struct('plate_dims', size(drugs_struct(1).layout), ...
    'treated_wells', true(size(drugs_struct(1).layout)), ...
    'Drugs', Drugs, 'Well_volume', drugs_struct(1).well_volume);

save(fullfile(pathstr, [name '_corr' ext]), 'Designs')
