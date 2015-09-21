% function Raw_data = Load_Plate_Data(plate_design)

function Raw_data = Load_Plate_Data(plate_design)

Raw_data = plate_design;

if isfield(plate_design,'FileSuffix')
    suffix = plate_design.FileSuffix;
else
    suffix = [];
end
fprintf('Loading: ');

for i = 1:plate_design.Nwells
    fprintf('%s ', plate_design.Wells{i,3});
    if mod(i,10)==0, fprintf('\n        ');end

    temp = load([plate_design.Folder ...
        plate_design.Wells{i,3} suffix '.mat']);
    if isfield(temp,'is_gap'), temp=rmfield(temp,'is_gap');end
    if strcmpi(Raw_data.RFPchannel, 'MOMP')
        if isfield(temp,'rfp')
            temp.momp = temp.rfp;
            temp = rmfield(temp,'rfp');
        end
    end
    if i>1
        Raw_data.rawdata(i) = orderfields(temp, Raw_data.rawdata(1));
    else
        Raw_data.rawdata(i) = temp;
    end

    if i==1
        Ntimepoints = size(Raw_data.rawdata(i).fret,1);
    end

    Raw_data.Ntrajectories(i) = size(Raw_data.rawdata(i).fret,2);
    assert_TrajLength(Raw_data.rawdata(i), plate_design.Wells{i,3}, ...
        [Ntimepoints Raw_data.Ntrajectories(i)]);



    Raw_data.WellLabel{i} = [];
    for j=find(Raw_data.Doses(i,:)>0)
        if ~isempty(Raw_data.WellLabel{i})
            Raw_data.WellLabel{i} = [Raw_data.WellLabel{i} ' + '];
        end
        Raw_data.WellLabel{i} = [ Raw_data.WellLabel{i} Raw_data.Drugs{j}];
        if ~strcmpi(Raw_data.Units{j},'inf')
            Raw_data.WellLabel{i} = [ Raw_data.WellLabel{i} ...
                ' ' num2str(Raw_data.Doses(i,j)) Raw_data.Units{j}];
        end
    end
end
Raw_data.Ntimepoints = Ntimepoints;
fprintf('\ndone\n\n');

end

function assert_TrajLength(rawdata, well, sizes)

for field = setdiff(fields(rawdata)',{'frmidx' 'dx' 'dy' 'rfp0'})
    assert(all(size(rawdata.(field{:}))==sizes), ...
        'Not-consistent matrix size for well %s and field %s',...
        well, field{:})
end
end
