function Table_PlateLayout_toXLSX(t_layout, filename)
% ExportMap(t_layout, filename)
%
%   label       treatment label
%   layout      16x24 matrix with plate concentration (in uM)
%

if ~exist('filename','var')
    filename = 'temp.xlsx';
end

% find the different plates
if ismember('Plate',t_layout.Properties.VariableNames)
    PlateNames = unique(t_layout.Plate);
    nPlate = length(PlateNames);
else
    PlateNames = [];
    nPlate = 1;
end

% find the well identity
if ismember('Well', t_layout.Properties.VariableNames)
    if ~iscategorical(t_layout.Well)
        t_layout.Well = categorical(t_layout.Well);
    end
elseif ismember('Row', t_layout.Properties.VariableNames) && ...
        ismember('Column', t_layout.Properties.VariableNames)
    Well = horzcat(cellstr(t_layout.Row), cellfun2(@(x) num2str(x,'%02i'), ...
        num2cell(t_layout.Column)));
    Well = categorical(Well);
    t_layout = [t_layout table(Well)];
else
    error('Missing Well or column/Row labels')
end

% for iCol = setdiff(t_layout.Properties.VariableNames, {'Well' 'Column' 'Row' 'Plate'})
%     if iscategorical(t_layout.(iCol{:}))
%         t_layout.(iCol{:}) = cellstr(t_layout.(iCol{:}));
%     end
% end


for iP = 1:nPlate

    output = cell(17,25);

    if ~isempty(PlateNames)
        if iscellstr(PlateNames)
            t_plate = t_layout(strcmp(t_layout.Plate,PlateNames{iP}),:);
            output{1,1} = PlateNames{iP};
        else
            t_plate = t_layout(t_layout.Plate==PlateNames(iP),:);

            output{1,1} = sprintf('Plate %i', PlateNames(iP));
        end
    else
        output{1,1} = 'Plate layout';
        t_plate = t_layout;
    end


    letter_label = unique(t_plate.Row,'sorted');
    num_label = unique(t_plate.Column,'sorted');


    output(2:(1+length(letter_label)),1) = letter_label;
    output(1,2:(1+length(num_label))) = num2cell(num_label');


    for iL = 1:length(letter_label)
        for iN = 1:length(num_label)
            well = [letter_label{iL} num2str(num_label(iN),'%02i')];
            temp = t_plate(t_plate.Well==well, setdiff(t_plate.Properties.VariableNames, ...
                {'Well' 'Column' 'Row'}));

            str = '';
            for i=1:length(temp.Properties.VariableNames)
                if ~isempty(temp.(temp.Properties.VariableNames{i})) && ...
                        ~isequal(temp.(temp.Properties.VariableNames{i}),0)
                    if ~isempty(str)
                        str = sprintf('%s\n', str);
                    end
                    str = sprintf('%s%s=%s',str,temp.Properties.VariableNames{i},...
                        ConvertToStr(temp.(temp.Properties.VariableNames{i})));
                end
            end


            output{1+iL, 1+iN} = str;

        end
    end
    sheet = output{1,1};

    if exist(filename, 'file')
        [~,sheets] = xlsfinfo(filename);
        if ismember(sheet, sheets)
            fprintf('!! Deleting sheet %s in file %s\n', sheet, filename)
            [~,~,dumb] = xlsread(filename, sheet);
            xlswrite(filename, cell(size(dumb)), sheet);
        end
    end

    xlswrite(filename, output, sheet)


end
