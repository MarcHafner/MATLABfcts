function [DrugNames, HMSLids] = splitHMSLid(DrugNames)

iscat = iscategorical(DrugNames);
if iscat
    DrugNames = cellstr(DrugNames);
    DrugNames(strcmp(DrugNames, '<undefined>')) = {''};
end

HMSLids = cellcell2cellstr( ...
    cellfun2(@(x) regexp(x, '(HMSL\d*)', 'tokens'), DrugNames));

%%
DrugNames = cellcell2cellstr(cellfun2(@(x) regexp(x, 'HMSL\d*', 'split'), DrugNames));
DrugNames = cellfun2(@strcat,DrugNames);
DrugNames = cellcell2cellstr(cellfun2(@(x) regexp(x, '[\s,]$', 'split'), DrugNames));
DrugNames = cellfun2(@strcat,DrugNames);

if iscat
    HMSLids = categorical(HMSLids);
    HMSLids(isundefined(HMSLids)) = '-';
    DrugNames = categorical(DrugNames);
    DrugNames(isundefined(DrugNames)) = '-';
end
