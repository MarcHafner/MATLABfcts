
function matrix2tsv(matrix, filename, labels, format)
% MatrixTotsv(matrix, filename, labels, format)
%
%       matrix:     2D matrix of numeric values
%       filename:   name of the export file
%       labels:     cellstr, or 2x1 cell of cellstr for dimensions 1 and 2
%       format:     format for writing numeric values (default: %.4g)
%

if ~exist('format','var')
    format = '%.4g';
end

f = fopen(filename,'w');
if exist('labels','var')
    if iscellstr(labels)
        if length(labels)==size(matrix,1)
            if length(labels)==size(matrix,2)
                warnprintf('Assuming labels to match first dimension')
            end
            labels = {labels {}};
        elseif length(labels)==size(matrix,2)
            labels = {{} labels};
        else
            error('labels has not the right length')
        end
    elseif iscell(labels)
        if length(labels)==2
            assert(all(size(matrix)==cellfun(@length, labels)), ...
                'wrong size for the labels (%i, %i), matrix is (%i, %i)', ...
                cellfun(@length, labels), size(matrix))
            for i=1:2
                if ~iscellstr(labels{i}) && isvector(labels{i})
                    if iscategorical(labels{i})
                        labels{i} = cellstr(labels{i});
                    else
                        labels{i} = num2cellstr(labels{i}, format);
                    end
                end
            end
        else
            error('labels has not the right format')
        end
    elseif isvector(labels)
        if length(labels)==size(matrix,1)
            if length(labels)==size(matrix,2)
                warnprintf('Assuming labels to match first dimension')
            end
            labels = {num2cellstr(labels, format) {}};
        elseif length(labels)==size(matrix,2)
            labels = {{} num2cellstr(labels, format)};
        else
            error('labels has not the right length')
        end
    else
        error('labels has not the right format')
    end
else
    labels = { {} {} };
end

if ~isempty(labels{2})
    if ~isempty(labels{1})
        fprintf(f, '\t');
    end
    fprintf(f, strjoin(labels{2},'\t'));
    fprintf(f, '\n');
end

for i=1:size(matrix,1)
    if ~isempty(labels{1})
        fprintf(f, [labels{1}{i} '\t']);
    end
    fprintf(f, format, matrix(i,1));
    fprintf(f, ['\t' format], matrix(i,2:end));
    fprintf(f, '\n');
end

pause(.2)
fclose(f);
