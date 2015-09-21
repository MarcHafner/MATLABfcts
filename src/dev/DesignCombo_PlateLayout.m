function drugs_struct = DesignCombo_PlateLayout(PrimaryDrugs, PrimaryDoses, PrimarySingleDoses, ...
    SecondaryDrugs, SecondaryDoses, SecondarySingleDoses, randomize_seed, edge_ctrl, nominal_conc)
% function drugs_struct = plate_design_combo(PrimaryDrugs, PrimaryDoses, PrimarySingleDoses, ...
%    SecondaryDrugs, SecondaryDoses, SecondarySingleDoses, randomize_seed, edge_ctrl, nominal_conc)
%
%
%   Doses are in uM
%   nominal_conc in mM (stock)
%   default volume is 60 ul
%   minimal dispensed volume is 20 pl, step of 10pl
%   maximal dispensed volume is 120 nl to have a dilution of 500-fold minimum

if ~exist('edge_ctrl','var') || isempty(edge_ctrl)
    edge_ctrl = true;
end

if randomize_seed>1e3
    s = RandStream('mt19937ar','Seed',mod(randomize_seed,1e9));
    RandStream.setGlobalStream(s);
end

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

drugs_struct = struct('name', [PrimaryDrugs;SecondaryDrugs], ...
    'nominal_conc', nominal_conc, ...
    'layout',zeros(16,24), 'well_volume', well_volume);


total_cnt = length(PrimaryDrugs)*length(PrimarySingleDoses) + ...
    length(SecondaryDrugs)*length(SecondarySingleDoses) + ...
    numel(PrimaryDoses)*numel(SecondaryDoses);


for iPD = 1:length(PrimaryDrugs)
    drugs_struct(iPD).Primary = true;
    drugs_struct(iPD).Doses = round_Conc(PrimaryDoses(iPD,:), ...
        drugs_struct(iPD).nominal_conc);
    drugs_struct(iPD).SingleDoses = round_Conc(PrimarySingleDoses(iPD,:), ...
        drugs_struct(iPD).nominal_conc);
    assert(all(ismember(PrimaryDoses(iPD,:), PrimarySingleDoses(iPD,:))),...
        'Dose in the combo not found in the titration for %s',PrimaryDrugs{iPD})
end
for iSD = 1:length(SecondaryDrugs)
    drugs_struct(length(PrimaryDrugs)+iSD).Primary = false;
    drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses = ...
        round_Conc(SecondarySingleDoses(iSD,:),drugs_struct(length(PrimaryDrugs)+iSD).nominal_conc);
    drugs_struct(length(PrimaryDrugs)+iSD).Doses = ...
        round_Conc(SecondaryDoses(iSD,:),drugs_struct(length(PrimaryDrugs)+iSD).nominal_conc);
    if size(SecondarySingleDoses,2)>0
        assert(all(ismember(SecondaryDoses(iSD,:), SecondarySingleDoses(iSD,:))),...
            'Dose in the combo not found in the titration for %s',SecondaryDrugs{iSD})
    end
end


fprintf('Number of test wells: %i\n',total_cnt);
assert(total_cnt<384, 'Too many well used (%i)', total_cnt)

ctrl_cnt = 384 - total_cnt;

if edge_ctrl
    assert(ctrl_cnt>=12, 'For controls on the edge, need at least 12 controls')
    ctrlpos = zeros(16,24);

    if ctrl_cnt<14 % only 6 controls on the edge
        disp('Only 6 control on the edge, would be better with at least 14 controls')
        ctrlpos([1 4 12 16],[7 18])=1;
        ctrlpos(8,[1 7 18 24])=1;
    else % 2 controls on each edge, 6 regularly spread in the middle.
        ctrlpos([1 4 12 16],[7 18])=1;
        ctrlpos([5 11],[1 12 24])=1;
    end

    if ctrl_cnt>=20 % put the corners as control (should be then discarded)
        ctrlpos([1 end], [1 end]) = 1;
    end


    if ctrl_cnt>sum(ctrlpos(:))
        % complete the number of control with randomized positions
        temp = reshape(randperm(16*24),16,24);
        temp(ctrlpos==1) = 384;
        cutoff = sort(temp(:));
        cutoff = cutoff(ctrl_cnt-sum(ctrlpos(:)));
        ctrlpos(temp<=cutoff) = 1;
    end

    ctrlidx = find(ctrlpos);
    assert(ctrl_cnt==length(ctrlidx));
else
    ctrlidx = find(reshape(randperm(16*24),16,24)<=(ctrl_cnt));
end

if ~exist('randomize','var') || randomize_seed>0
    trtidx = randperm(384);
    trtidx = trtidx(~ismember(trtidx,ctrlidx));
else
    trtidx = setdiff(1:384,ctrlidx);
end

letter_label = 'A':'P';

fprintf('Control wells (%i):\n', ctrl_cnt);
for i=1:length(ctrlidx)
    fprintf('\t%s%i', letter_label(mod(ctrlidx(i)-1,16)+1), ceil(ctrlidx(i)/16));
end
disp(' ')



cnt = 1;


for iPD = 1:length(PrimaryDrugs)

    % titration
    for iPSDo = 1:length(drugs_struct(iPD).SingleDoses)
        drugs_struct(iPD).layout(trtidx(cnt)) = drugs_struct(iPD).SingleDoses(iPSDo);
        cnt = cnt+1;
    end

    for iSD = 1:length(SecondaryDrugs)

        % combos
        for iPDo = 1:length(drugs_struct(iPD).Doses)
            for iSDo = 1:length(drugs_struct(length(PrimaryDrugs)+iSD).Doses);



                drugs_struct(iPD).layout(trtidx(cnt)) = drugs_struct(iPD).Doses(iPDo);

                drugs_struct(length(PrimaryDrugs)+iSD).layout(trtidx(cnt)) = ...
                    drugs_struct(length(PrimaryDrugs)+iSD).Doses(iSDo);
                vol = well_volume*((drugs_struct(iPD).Doses(iPDo)/drugs_struct(iPD).nominal_conc)+...
                    (drugs_struct(length(PrimaryDrugs)+iSD).Doses(iSDo)/...
                    drugs_struct(length(PrimaryDrugs)+iSD).nominal_conc))/1e3;
                if vol>max_volumedrop
                    warning('Too large volume used (%.2f%%) for combo %s(%.2fuM) - %s(%.2fuM)', ...
                        100*vol/well_volume, drugs_struct(iPD).name, drugs_struct(iPD).Doses(iPDo), ...
                        drugs_struct(length(PrimaryDrugs)+iSD).name, ...
                        drugs_struct(length(PrimaryDrugs)+iSD).Doses(iSDo))
                end

                cnt = cnt+1;
            end
        end

    end
end


for iSD = 1:length(SecondaryDrugs)
    % titration
    for iSSDo = 1:length(drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses)
        drugs_struct(length(PrimaryDrugs)+iSD).layout(trtidx(cnt)) = ...
            drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses(iSSDo);
        cnt = cnt+1;
    end

end



assert( cnt==length(trtidx)+1)

for iPD = 1:length(PrimaryDrugs)
    assert(sum(drugs_struct(iPD).layout(:)~=0) == (length(drugs_struct(iPD).Doses)*...
        sum((~[drugs_struct.Primary]).*cellfun(@length, {drugs_struct.Doses}))+...
        length(drugs_struct(iPD).SingleDoses)))
    % volume in nl from stock (conc in mM)
    drugs_struct(iPD).volume = 1e9*well_volume*sum(drugs_struct(iPD).layout(:))/...
        (1e3*drugs_struct(iPD).nominal_conc);
end

for iSD = 1:length(SecondaryDrugs)
    assert(sum(drugs_struct(length(PrimaryDrugs)+iSD).layout(:)~=0)==...
        (length(drugs_struct(length(PrimaryDrugs)+iSD).Doses)*...
        sum([drugs_struct.Primary].*cellfun(@length, {drugs_struct.Doses}))+...
        length(drugs_struct(length(PrimaryDrugs)+iSD).SingleDoses)))
    drugs_struct(length(PrimaryDrugs)+iSD).volume = ...
        1e9*well_volume*sum(drugs_struct(length(PrimaryDrugs)+iSD).layout(:))/...
        (1e3*drugs_struct(length(PrimaryDrugs)+iSD).nominal_conc);
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
