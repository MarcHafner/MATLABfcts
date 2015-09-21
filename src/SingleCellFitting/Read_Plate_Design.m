% function plate_design = Read_Plate_Design(filename, tag, verbatim)
%
%   filename:   <path>/<filename>.tsv
%   [tag]:      tag for the experiment (default = filename)
%   [verbatim]: display read lines (default = false)
%
%   read the plate design in the filename file with format:
%
%             Folder	../data/20140306_Mapa_Fc/
%             Nwells	28
%             Experimentalist	Sam
%             Date	20140306
%             RFPchannel	MOMP
%             Timestep	5
%             Drugs	Mapa	Fc	TRAIL
%             Units	nM	-	ng/ml
%             Treatments
%             C04	0	0	0
%               .....
%             F10   10  1   0
%
%   NOTE: well name should correspond to the file names (C4 vs c4 vs C04).
%

function plate_design = Read_Plate_Design(filename, tag, verbatim)

if ~exist('tag','var') || isempty(tag)
    [~,tag] = fileparts(filename);
end

if ~exist('verbatim','var')
    verbatim = false;
end
[file, err] = fopen(filename,'r');
if file==-1
    error('%s: %s',err, filename);
end

plate_design.tag = tag;

while ~feof(file)
    line = fgetl(file);
    line = regexp(line,'\t','split');
    if verbatim,disp(['   ' cellstr2str(line,'~')]),end

    if ~isempty(strfind(line{1},'Treatm')) || ~isempty(strfind(line{1},'treatm'))
        disp('all input fields read')
        break,end
    idx = setdiff(find(~cellfun(@isempty,line)),1);
    if length(idx)==1
        plate_design.(line{1}) = line{2};
    else
        plate_design.(line{1}) = line(idx);
    end
end
assert(~feof(file),'Line ''Treatments'' not found!')
for field = {'Drugs' 'Nwells' 'Date' 'Timestep' 'RFPchannel'}
    assert(isfield(plate_design,field{:}), 'Need a ''%s'' field', field{:})
end

nDrugs = length(plate_design.Drugs);
if isfield(plate_design,'Units')
    assert(nDrugs==length(plate_design.Units), 'Need as many ''Units'' as ''Drugs''')
end

plate_design.Timestep = str2double(plate_design.Timestep);
plate_design.Nwells = str2double(plate_design.Nwells);
plate_design.Wells = cell(plate_design.Nwells,3);
plate_design.Doses = NaN(plate_design.Nwells,nDrugs);
Wellcnt = 0;

while ~feof(file)

    line = fgetl(file);
    line = regexp(line,'\t','split');
    if verbatim,disp(['   ' cellstr2str(line,'~')]),end

    Wellcnt = Wellcnt+1;
    plate_design.Wells(Wellcnt,:) = {line{1}(1) str2double(line{1}(2:end)) line{1}};
    plate_design.Doses(Wellcnt,:) = cell2mat(cellfun2(@str2double,line(2:end)));

    if feof(file)
        assert(Wellcnt == plate_design.Nwells, ...
            'Number of treatments (%i) inconsistent with specified number of wells (%i)', ...
            Wellcnt, plate_design.Nwells)
        disp('all treatments read')
        break
    end

end
fclose(file);

plate_design.axes = {unique(plate_design.Wells(:,1)) ...
    unique(cell2mat(plate_design.Wells(:,2)))};

plate_design.DrugConc = NaN(length(plate_design.axes{1}), ...
    length(plate_design.axes{2}), nDrugs);

for i=1:length(plate_design.axes{1})
    for j=1:length(plate_design.axes{2})
        idx = find( strcmp(plate_design.Wells(:,1),plate_design.axes{1}{i}) & ...
            cell2mat(plate_design.Wells(:,2))==plate_design.axes{2}(j));
        assert(length(idx)<2);
        if ~isempty(idx)
            plate_design.DrugConc(i,j,:) = plate_design.Doses(idx,:);
        end
    end
end

figure(991);clf
for i=1:nDrugs
    subplot(floor(sqrt(nDrugs)), ceil(nDrugs/floor(sqrt(nDrugs))), i)
    imagesc(plate_design.DrugConc(:,:,i))
    colormap(jet(256))
    colorbar;
    title([plate_design.Drugs{i} ' ' plate_design.Units{i}])
    set(gca,'xtick',1:length(plate_design.axes{2}),'xticklabel', ...
        plate_design.axes{2}, 'ytick',1:length(plate_design.axes{1}),'yticklabel', ...
        plate_design.axes{1})
end
