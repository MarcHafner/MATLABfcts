function lineArray = read_mixed_tsv(fileName,delimiter)

if ~exist('delimiter','var')
    delimiter = '\t';
end

fid = fopen(fileName,'r');      % Open the file
lineArray = cell(100,1);        % Preallocate a cell array (ideally slightly
                                %   larger than is needed)
lineIndex = 1;
nextLine = fgetl(fid);
while ~isequal(nextLine,-1)         % Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;
    lineIndex = lineIndex+1;
    nextLine = fgetl(fid);
end
fclose(fid);                    % Close the file

lineArray = lineArray(1:lineIndex-1);               % Remove empty cells, if needed
for iLine = 1:lineIndex-1
    lineData = textscan(lineArray{iLine},'%s',...   % Read strings
        'Delimiter',delimiter);
    lineData = lineData{1};
    if strcmp(lineArray{iLine}(end),delimiter)      % Account for when the line
        lineData{end+1} = '';                       %   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  % Overwrite line data
end
end
