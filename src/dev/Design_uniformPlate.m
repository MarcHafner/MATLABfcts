function drugs_struct = Design_uniformPlate(Drugs, Doses,  nominal_conc)
% function drugs_struct = plate_design_combo(PrimaryDrugs, PrimaryDoses, PrimarySingleDoses, ...
%    SecondaryDrugs, SecondaryDoses, SecondarySingleDoses, randomize_seed, edge_ctrl, nominal_conc)
%
%
%   Doses are in uM
%   nominal_conc in mM (stock)
%   default volume is 60 ul
%   minimal dispensed volume is 20 pl, step of 10pl
%   maximal dispensed volume is 120 nl to have a dilution of 500-fold minimum


% well volume
well_volume = 6e-5;
min_volumedrop = 1e-11; % minimum step of the drop is 10pl for the
max_volumedrop = well_volume/500; % need a dilution of 500-fold minimum

% stock conc in mM
if ~exist('nominal_conc','var')
    nominal_conc = 10;
elseif isvector(nominal_conc)
    nominal_conc = num2cell(ToColumn(nominal_conc));
end

% conc in mM

drugs_struct = struct('name', Drugs, ...
    'nominal_conc', nominal_conc, ...
    'layout',zeros(16,24), 'well_volume', well_volume);

for iD = 1:length(Drugs)
    drugs_struct(iD).SingleDoses = round_Conc(Doses(iD), ...
        drugs_struct(iD).nominal_conc);

    %%% a specific layout could be set as an input ; uniform for now.
    drugs_struct(iD).layout(:) = drugs_struct(iD).SingleDoses;

    % volume in nl from stock (conc in mM)
    drugs_struct(iD).volume = 1e9*well_volume*sum(drugs_struct(iD).layout(:))/...
        (1e3*drugs_struct(iD).nominal_conc);
end


    function new_Doses = round_Conc(old_Doses, nominal_conc)

        dose_step = 1e3*nominal_conc *min_volumedrop/well_volume;
        dose_max = 1e3*nominal_conc *max_volumedrop/well_volume;
        temp_d = old_Doses;
        if any(temp_d<2*dose_step)
            disp(['!! Doses below minimal conc (' num2str(2*dose_step) ')'])
            disp(['       at doses : ' num2str(temp_d(temp_d<2*dose_step))])
            temp_d = max(temp_d, 2*dose_step);
        end
        if any(temp_d<50*dose_step)
            idx = temp_d<50*dose_step;
            disp(['former doses: ' num2str(temp_d(idx))])
            temp_d(idx) = dose_step*round(temp_d(idx)/dose_step);
            disp(['      now as: ' num2str(temp_d(idx))])
        end
        if any(temp_d>dose_max)
            disp(['!! Doses above maximal conc (' num2str(dose_max) ')'])
            disp(['       at doses : ' num2str(temp_d(temp_d>dose_max))])
            temp_d = min(temp_d, dose_max);
        end
        new_Doses = temp_d;
    end

end
