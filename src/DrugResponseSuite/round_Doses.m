
function new_Doses = round_Doses(old_Doses, nominal_conc, DName, tag, ...
    min_volume, step_volume, max_volume, well_volume)
% new_Doses = round_Doses(old_Doses, nominal_conc, DName, tag, ...
%    min_volume, step_volume, max_volume, well_volume)
% adjust doses to match the specifications of the D300
%

min_dose = nominal_conc *min_volume/well_volume;
step_dose = nominal_conc *step_volume/well_volume;
max_dose = nominal_conc *max_volume/well_volume;
temp_Doses = old_Doses(old_Doses>0);
old_NZ_Doses = temp_Doses;
new_Doses = old_Doses;
fprintf('Processing %s (%s doses):\n', DName, tag)
if any(temp_Doses<min_dose)
    warnprintf('%i dose(s) below minimal conc (set to %.2g)', ...
        sum(temp_Doses<min_dose), min_dose);
    temp_Doses = max(temp_Doses, min_dose);
end
if any(temp_Doses>max_dose)
    warnprintf('%i dose(s) above max conc (set to %.2g)', ...
        sum(temp_Doses>max_dose), max_dose);
    temp_Doses = min(temp_Doses, max_dose);
end
if any( ((mod(temp_Doses,step_dose)./temp_Doses)>.02) & (temp_Doses>1.02*min_dose) )
    idx = ((mod(temp_Doses,step_dose)./temp_Doses)>.02) & (temp_Doses>1.02*min_dose);
    fprintf('\tNote: %i doses round down because more than 2%% difference\n', ...
        sum(idx));
    % case of the lowest doses
    idx = (temp_Doses<mean([min_dose, step_dose])) & (temp_Doses>1.02*min_dose);
    temp_Doses(idx) = min_dose;
    idx = (temp_Doses>mean([min_dose, step_dose])) & (temp_Doses<step_dose/1.02);
    temp_Doses(idx) = step_dose;    
    % other cases
    idx = ((mod(temp_Doses,step_dose)./temp_Doses)>.02) & (temp_Doses>1.02*step_dose);
    temp_Doses(idx) = step_dose*round(temp_Doses(idx)/step_dose);
end
if any(temp_Doses~=old_NZ_Doses)
    fprintf('\tChanged doses:')
    for iDoses = ToRow(find(temp_Doses~=old_NZ_Doses))
        fprintf('\t%-7.2g', old_NZ_Doses(iDoses))
    end
    fprintf('\n\t    (%3i)  -->', sum(temp_Doses~=old_NZ_Doses))
    for iDoses = ToRow(find(temp_Doses~=old_NZ_Doses))
        fprintf('\t%-7.2g', temp_Doses(iDoses))
    end
    fprintf('\n')
end

new_Doses(old_Doses>0) = temp_Doses;
