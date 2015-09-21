function compiled_data = compile_ExpData(folder, SelectedDrug, coDrugs, savefile)
% compiled_data = compile_ExpData(folder, SelectedDrug, coDrug, savefile)
%
%   Example: compile all TRAIL data with ABT or Bortez
%     folder = './datasets/';
%     SelectedDrug = 'TRAIL';
%     coDrugs = {'Bortezomib' 'ABT'};
%     savefile = 'allTRAIL.mat';
%

saving = true;
if ~exist('savefile','var')
    saving = false;
elseif isnumeric(savefile) || islogical(savefile)
    savefile = ['all' SelectedDrug '.mat'];
end

%%
filelist = dir([folder '*.mat']);
filelist = {filelist.name};

Selectedfields = [ ...
    SelectedDrug ReplaceName(coDrugs) ...
    {'Experimentalist' 'Date' 'Legend'}];
Fitfields = {'r2' 'k' 'Delay' 'badidx'};
Fatefields = {'Surviving' 'MompTime' 'FRETSurviving' 'FRETMompTime' 'FRETPreMompTime' ...
    'preMompMaxActivity' 'preMompMaxActivityTime' 'MaxIC' 'MaxActivity' ...
    'MaxActivityTime' 'MaxFRETPreFRETMomp' 'FRETMompTimeIdx' 'FRETPreMompTimeIdx'};
Trajfields = {'ICTraj' 'C8Activity'};

Rawfields = {'fret' 'area' 'edge'};
rawTrajfields = {'Prob'};

temp = cell(1,length(Selectedfields));
warning('off', 'MATLAB:codetools:ModifiedVarnames')
warning('off','MATLAB:table:RowsAddedExistingVars');

compiled_data = struct('ExpKey',table( temp{:},'VariableNames', Selectedfields),...
    'fits',struct,'Traj',struct,'Tracking',table,'Fate',struct,'rawdata',struct);

cnt = 0;
allrfp0 = cell(0,1);
for iF=1:length(filelist)
    fprintf('Loading file %s: ', filelist{iF});

    load([folder filelist{iF}])

    Drugidx = find(strcmp(SelectedDrug, Fitted_data.Drugs));
    if isempty(Drugidx)
        fprintf(' no %s treatment -> skip\n\n', SelectedDrug);
        continue
    end
    coDrugidx = find(ismember(Fitted_data.Drugs, coDrugs));
    notDrugidx = setdiff(1:length(Fitted_data.Drugs), [Drugidx, coDrugidx]);
    notDrugs = Fitted_data.Drugs(notDrugidx);

    if isempty(notDrugidx)
        Widx = find(Fitted_data.Doses(:,Drugidx)~=0);
        notDrugs = {''};
    else
        Widx = find(Fitted_data.Doses(:,Drugidx)~=0 & ...
            all(Fitted_data.Doses(:,notDrugidx)==0,2));
    end
    if isempty(Widx)
        fprintf(' no %s treatment without another drug -> skip\n\n', SelectedDrug);
        continue
    end

    fprintf('\n  Drug used: %s\n', cellstr2str(Fitted_data.Drugs))
    fprintf('\n  found %i wells with %s (and ev. %s):\n', length(Widx), SelectedDrug, ...
        cellstr2str(coDrugs));

    for iW=Widx'
        cnt = cnt+1;

        doses = Fitted_data.Doses(iW,Drugidx);
        for i=1:length(coDrugs)
            idx = find(strcmp(Fitted_data.Drugs, coDrugs{i}));
            if isempty(idx)
                doses(i+1) = 0;
            else
                doses(i+1) =  Fitted_data.Doses(iW,idx);
            end
        end
        fprintf('\t%s=%.1f, [%s] = [%s] ; [%s] = [%s]\n', SelectedDrug, Fitted_data.Doses(iW,Drugidx), ...
            cellstr2str(coDrugs), num2str(doses(2:end)), cellstr2str(notDrugs),...
            num2str(Fitted_data.Doses(iW,notDrugidx)));

        compiled_data.ExpKey = [ compiled_data.ExpKey
            cell2table([num2cell(doses) ...
            Fitted_data.Experimentalist Fitted_data.Date Fitted_data.WellLabel{iW}], ...
            'VariableNames', Selectedfields)];

        for ifield = Fitfields
            compiled_data.fits(cnt).(ifield{:}) = ...
                Fitted_data.fit(iW).(ifield{:});
        end

        compiled_data.Tracking(cnt,'Tracking') = {Fitted_data.Ntrajectories(iW)};
        for iK=Fitted_data.Legends.ErrorType.keys
            VarName = internal.matlab.codetools.genvalidnames({Fitted_data.Legends.ErrorType(iK{:})});
            compiled_data.Tracking(cnt,VarName) = ...
                {sum(Fitted_data.Traj(iW).ErrorType==iK{:})};
        end

        compiled_data.Fate(cnt).('NGoodTraj') = Fitted_data.NgoodTraj(iW);
        for ifield = Fatefields
            if islogical(Fitted_data.stats(iW).(ifield{:}))
                compiled_data.Fate(cnt).(['N' ifield{:}]) = ...
                    sum(Fitted_data.stats(iW).(ifield{:}));
            end
            compiled_data.Fate(cnt).(ifield{:}) = ...
                Fitted_data.stats(iW).(ifield{:});
        end

        for ifield = Trajfields
            compiled_data.Traj(cnt).(ifield{:}) = ...
                Fitted_data.goodTraj(iW).(ifield{:});
        end
        compiled_data.Traj(cnt).T = Fitted_data.T;


        for ifield = Rawfields
            compiled_data.rawdata(cnt).(ifield{:}) = ...
                Fitted_data.rawdata(iW).(ifield{:});
        end

        for ifield = rawTrajfields
            compiled_data.rawdata(cnt).(ifield{:}) = ...
                Fitted_data.Traj(iW).(ifield{:});
        end
        compiled_data.rawdata(cnt).SelectedCells = ...
            Fitted_data.Traj(iW).goodtrack;


        if isfield(Fitted_data,'rfp0')
            allrfp0{cnt} = Fitted_data.rfp0{iW};
        else
            allrfp0{cnt} = {};
        end

    end
    fprintf('\n');
end
if any(~cell2mat(cellfun2(@isempty,allrfp0)))
    compiled_data.rfp0 = allrfp0;
end

warning('on', 'MATLAB:codetools:ModifiedVarnames')
warning('on', 'MATLAB:table:RowsAddedExistingVars');

save(savefile, 'compiled_data')
