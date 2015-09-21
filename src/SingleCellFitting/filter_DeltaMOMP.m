function filtered_data = filter_DeltaMOMP(compiled_data, cutoffs)

if ~exist('cutoffs','var')
    cutoffs = [-90 90];
end
get_newfigure(9999)
filtered_data = compiled_data;
lD = height(compiled_data.ExpKey);
for iE = 1:lD

    subplot(floor(sqrt(lD)),ceil(lD/floor(sqrt(lD))),iE)

    Goodidx = compiled_data.Fate(iE).Surviving | ...
        ( ((compiled_data.Fate(iE).MompTime-compiled_data.Fate(iE).FRETMompTime)<cutoffs(2) & ...
        (compiled_data.Fate(iE).MompTime-compiled_data.Fate(iE).FRETMompTime)>cutoffs(1))  & ...
        (compiled_data.Fate(iE).MompTime>100 | compiled_data.fits(iE).k>0)) ;

    hist(compiled_data.Fate(iE).MompTime-compiled_data.Fate(iE).FRETMompTime,...
        (cutoffs(1)-40):10:(cutoffs(2)+40));
    hold on
    plot(cutoffs([1 1]), ylim, '-r')
    plot(cutoffs([2 2]), ylim, '-r')
    xlim([(cutoffs(1)-50) (cutoffs(2)+50)])

    nCells = length(compiled_data.fits(iE).r2);
    fprintf('%i (%.1f%%) cells selected from %i (%s ; %s)\n', sum(Goodidx), ...
        100*(sum(Goodidx)/length(Goodidx)), length(Goodidx), ...
        compiled_data.ExpKey.Legend{iE}, compiled_data.ExpKey.Date{iE});
    if (sum(Goodidx)/length(Goodidx))<.8
        warning('Less than 80%% of cells selected for %s',...
            compiled_data.ExpKey.Legend{iE})
    end

    filtered_data.MOMPfilter{iE} = Goodidx;

    for f = fieldnames(compiled_data)'
        if isstruct(filtered_data.(f{:}))
        for f2 = fieldnames(filtered_data.(f{:})(iE))'
            idx = find(size(compiled_data.(f{:})(iE).(f2{:}))==nCells);
            assert(length(idx)<2)
            if isempty(idx),continue,end
            if idx==1
                filtered_data.(f{:})(iE).(f2{:}) = ...
                    compiled_data.(f{:})(iE).(f2{:})(Goodidx,:);
            elseif idx==2
                filtered_data.(f{:})(iE).(f2{:}) = ...
                    compiled_data.(f{:})(iE).(f2{:})(:,Goodidx);
            else
                error('Mismatched dimensions for field %s', f{:});
            end
        end
        elseif iscell(filtered_data.(f{:})) && ...
                length(filtered_data.(f{:}))==height(filtered_data.ExpKey) &&...
                ~isempty(filtered_data.(f{:}){iE})
            filtered_data.(f{:}){iE} = filtered_data.(f{:}){iE}(Goodidx);
        end
    end

end
