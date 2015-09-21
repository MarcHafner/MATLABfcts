function ExportMap_D300(drugs_struct, filename, sheet)
% ExportMap(drug_struct, filename, sheet)
%
% drug_struct.
%   name        drug name
%   conc        cateridge concentration (in mM)
%   layout      16x24 matrix with plate concentration (in uM)
%

maxDrugName = 10; % length of drug name for layout sheet

%output = cell(19*length(drugs_struct),25);
total_output = cell(18+4*length(drugs_struct),25);
letter_label = 'A':'P';
num_label = 1:24;

total_output(2:17,1) = num2cell(letter_label);
total_output(1,2:end) = num2cell(num_label');

for i=1:length(drugs_struct)
%     temp_out = cell(19,25);
%     temp_out{1,1} = drugs_struct(i).name;
%     temp_out{1,2} = num2str(drugs_struct(i).nominal_conc);
%     temp_out{1,3} = 'mM';
%
%     temp_out(3:end-1,1) = num2cell(letter_label);
%     temp_out(2,2:end) = num2cell(num_label');
%     temp_out(3:end-1,2:end) = num2cell(drugs_struct(i).layout);
%
%     output((1:19)+(i-1)*19,:) = temp_out;

    for i1 = 1:16
        for i2 = 1:24
            if drugs_struct(i).layout(i1,i2)>0
                if isempty(total_output{i1+1,i2+1})
                    total_output{i1+1,i2+1} = sprintf('%s, %.2fuM', ...
                        drugs_struct(i).name(1:min(end,maxDrugName)), ...
                        drugs_struct(i).layout(i1,i2));
                else
                    total_output{i1+1,i2+1} = sprintf('%s\n%s, %.2fuM', ...
                        total_output{i1+1,i2+1}, drugs_struct(i).name(1:min(end,maxDrugName)),...
                        drugs_struct(i).layout(i1,i2));
                end
            end
        end
    end

    temp = {drugs_struct(i).name 'Stock (mM)=' drugs_struct(i).nominal_conc ...
        'Volume (nl)=' drugs_struct(i).volume 'Well volume (ul)=' drugs_struct(i).well_volume*1e6 };
    if isfield(drugs_struct(i),'Doses')
        temp(2,1:(length(drugs_struct(i).Doses)+1)) = [{'Combo Doses (uM)='} ...
            num2cell(drugs_struct(i).Doses)];
    end
    if isfield(drugs_struct(i),'SingleDoses')
        temp(3,1:(length(drugs_struct(i).SingleDoses)+1)) = [{'Single Doses (uM)='} ...
            num2cell(drugs_struct(i).SingleDoses)];
    end
    total_output(18+(i-1)*4+(1:3),1:size(temp,2)) = temp;

end


if isfield(drugs_struct,'seeding')
    for i1 = 1:16
        for i2 = 1:24
            total_output{i1+1,i2+1} = sprintf('%s\nseed=%.0f', ...
                total_output{i1+1,i2+1},drugs_struct(1).seeding(i1,i2));
        end
    end
end

if isfield(drugs_struct,'Ligand')
    for i1 = 1:16
        for i2 = 1:24
            total_output{i1+1,i2+1} = sprintf('%s\n%s, %.2f', ...
               total_output{i1+1,i2+1}, drugs_struct(1).Ligand, drugs_struct(1).Ligand_layout(i1,i2));
        end
    end
end


if ~exist('filename','var')
    filename = 'temp.tsv';
end
if ~exist('sheet','var')
    sheet = 'default';
end

tsvwrite(filename, total_output, sheet)
