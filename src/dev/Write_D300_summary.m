function Write_D300_summary(filename, experiment)

%%
Drugs = {};
Stock = [];
Volume = [];

output = ['Drug Name' 'HMSL id' 'Stock (mM)' ...
    strcat({experiment.name}, ' (ul)') 'Total volume (ul)'];

for iS = 1:length(experiment)

    for i=1:length(experiment(iS).plate_map)
        t_Drug = experiment(iS).plate_map(i).name;
        Didx = find(strcmp(Drugs, t_Drug));

        stock_conc = experiment(iS).plate_map(i).nominal_conc;
        volume = experiment(iS).plate_map(i).volume;
        volume_str = num2str(volume/1e3,'%.2f');
        if isempty(Didx)
            Drugs{end+1} = t_Drug;  %#ok<AGROW>
            Stock(end+1) = stock_conc;  %#ok<AGROW>
            Volume(end+1,iS) = volume;  %#ok<AGROW>
            if strfind(t_Drug, 'HMSL')
                parts = regexp(t_Drug,'HMSL','split');
                t_Drug = parts{1};
                drug_id = ['HMSL' parts{2}];
            else
                drug_id = '_';
            end
            output(end+1,1:3) = {t_Drug drug_id stock_conc};  %#ok<AGROW>
            output{end,iS+3} = volume_str;
        else
            assert(Stock(Didx)==stock_conc, ...
                'not the same stock conc for %s', Drugs{Didx});
            Volume(Didx,iS) = volume;  %#ok<AGROW>
            output{Didx+1,iS+3} = volume_str;
        end
    end
end

for i=1:length(Drugs)
    output{i+1,end} = num2str(ceil(sum(Volume(i,:)/100))/10,'%.1f');
end


%%
tsvwrite(filename, output, 'summary_drugs');
