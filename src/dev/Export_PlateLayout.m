function Export_PlateLayout(layout, label, filename, sheet)
% ExportMap(layout, label, filename, sheet)
%
%   label       treatment label
%   layout      16x24 matrix with plate concentration (in uM)
%


output = cell(19,25);
letter_label = 'A':'P';
num_label = 1:24;


output{1,1} = label;

output(3:end-1,1) = num2cell(letter_label);
output(2,2:end) = num2cell(num_label');
output(3:end-1,2:end) = num2cell(layout);


if ~exist('filename','var')
    filename = 'temp.xlsx';
end

if ~exist('sheet','var')
    if exist(filename, 'file')
        warning('Deleting file %s', filename)
        delete(filename)
    end

    xlswrite(filename,output,'layout')

else
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
