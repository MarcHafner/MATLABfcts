% example to process all the datasets from a set of plate design


folder = './plate_design/';

filelist = dir([folder '*.tsv']);
filelist = {filelist.name};

%%

report = cell(1,length(filelist));

parfor i=1:length(filelist)
    try
        output = LoadFitData([folder filelist{i}]);
        report{i} = [filelist{i} ' is done with tag: ' num2str(output)];
    catch err
        warning([filelist{i} ' has error: ' err.message])
        report{i} = [filelist{i} ' has error: ' err.message];
    end
end
