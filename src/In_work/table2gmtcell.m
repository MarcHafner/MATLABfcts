
function gmt_output = table2gmtcell(t_GeneSet, GeneSetColumn, GeneColumn, label)
% gmt_output = table2gmtcell(t_GeneSet, GeneSetColumn, GeneColumn, label)
%
%   Generate a cell array that has the format of the gmt tsv file. It can
%   be used as a input for the GSEA anlysis.
%
%   t_Families is a table with a column for Gene set/Family and a column
%   for the genes.
%
%   GeneSetColumn is to specify the name of the column for the set (default
%   is 'GeneSet' 'Family' 'GeneFamily')
%
%   GeneColumn is to specify the name of the column for the genes (default
%   is 'Gene' 'GeneName')
%
%   label is to specify a label (source) for the set
%
%

t_GeneSet = TableToCategorical(t_GeneSet);

if ~exist('GeneSetColumn','var') || isempty(GeneSetColumn)
    GeneSetColumn = intersect(varnames(t_GeneSet), ...
        {'Family' 'GeneFamily' 'GeneSet'});
end

if ~exist('GeneColumn','var') || isempty(GeneColumn)
    GeneColumn = intersect(varnames(t_GeneSet), {'Gene' 'GeneName'});
end

if ~exist('label','var') || isempty(label)
    label = 'User Input List';
end

%%
Sets = unique(t_GeneSet.(GeneSetColumn));

gmt_output = cell(length(Sets), 3);
for i=1:length(Sets)
    gmt_output{i,1} = char(Sets(i));
    gmt_output{i,2} = label;
    genes = t_GeneSet.(GeneColumn)(t_GeneSet.(GeneSetColumn)==Sets(i));
    gmt_output{i,3} = strjoin(cellstr(genes),'\t');
end
