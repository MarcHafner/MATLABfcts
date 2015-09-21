function t_out = splitHMSLid_table(t_in, varname)

if~exist('varname','var')
    varname = 'DrugName';
end

DrugNames = t_in.(varname);

[DrugNames, HMSLids] = splitHMSLid(DrugNames);


idx = find(strcmp(t_in.Properties.VariableNames, varname));

t_out = [t_in(:,1:(idx-1)) table(DrugNames, HMSLids, 'variableNames', ...
    {'DrugName', 'HMSLid'}) t_in(:,(idx+1):end)];
