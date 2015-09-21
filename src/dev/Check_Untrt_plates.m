timestamp = '20140220';

load(['L1000validation_results_' timestamp '.mat'])

CellLines = categorical(1,length(t_CL));
for i=1:length(t_CL)
    CellLines(i) = t_CL{i}.CellLine(1);
end

%%


for iCL=1:3
    get_newfigure(100+iCL,[100 100 700 550])

    t_untrt = t_CL{iCL}(t_CL{iCL}.Replicate<0,{'Row','Column','Nuclei_NumberOfObjects'});

    [Untrt_cnt, labels] = table_to_ndarray( t_untrt, 'Key', {'Row','Column'}, 'outer', 1);

    get_newaxes([.1 .1 .5 .75])
    imagesc(Untrt_cnt,[500 1500])
    colormap(Plotting_parameters.cmapWBr)
    colorbar('location','northoutside')

    edgeWells = t_untrt( t_untrt.Row==1 | t_untrt.Row==16 | ...
        t_untrt.Column==1 | t_untrt.Column==24, :);
    centerWells = setdiff(t_untrt,edgeWells);

    get_newaxes([.65 .1 .3 .5], 1)

    x = 100:20:1500;
    n = ksdensity(edgeWells.Nuclei_NumberOfObjects,x);
    h = plot(x,n);

    n = ksdensity(centerWells.Nuclei_NumberOfObjects,x);
    h(2) = plot(x,n,'r');

    legend(h, {'Edge' 'center'})

    title(sprintf('\\Delta=%.1f, p=%.3f', ...
        median(centerWells.Nuclei_NumberOfObjects)-median(edgeWells.Nuclei_NumberOfObjects), ...
        ranksum(centerWells.Nuclei_NumberOfObjects,edgeWells.Nuclei_NumberOfObjects)))


end
