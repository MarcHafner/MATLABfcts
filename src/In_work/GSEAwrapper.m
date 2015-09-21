
function [t_results, Genelist] = GSEAwrapper(Genelist, Geneset, varargin)
% [t_results, Genelist] = GSEAwrapper(Genelist, Geneset, varargin)
%
%   Inputs:
%   - Genelist:
%       - filename
%       - cell array or table (first column is gene symbol,
%               2nd column is weight)
%   - Geneset:
%       - tag (GObp, GOcc, GOmf, allGO, KEGG, Reactome, Biocarta,
%               CGN, CGP, OncoSig, BreastCGP)
%       - filename
%       - cell array formated as gmt: each rows has the format:
%           { 'Set name' 'source' 'geneA\tgeneB\t...'}
%   - options:
%       - Outputfolder (folder in which the GSEA output folder will be
%           saved; default is none, results are not saved)
%       - Nplot (number of sets to display for html output; needs an
%               Outputfolder)
%       - Outputname (renamed the GSEA default folder name; needs
%               'Outputfolder'; WARNING: override previous results!)
%       - label (for the GSEA output, default is geneset label)
%       - Randseed (value for randomization)
%       - scoring_scheme (as defined by GSEA: weighted, classic or weighted_p2)
%       - set_min (minimal number of genes in a set; default=15)
%       - set_max (maximal number of genes in a set; default=500)
%       - GSEAfolder (where the GSEA.jar file is stored;
%               default is WorkFolder/GSEA_java)
%
%
%   Output:
%   - t_results: table with results for the selected sets
%   - Genelist : ordered gene lists
%
% requires the following files from GSEA to be stored in a unique folder:
%   - gsea2-2.1.0.jar
%   - gene set list such as (default), that can be called with [xx]
%       - c5.mf.v4.0.symbols.gmt [GOmf]
%       - c5.cc.v4.0.symbols.gmt [GOcc]
%       - c5.bf.v4.0.symbols.gmt [GObf]
%       - c2.cp.kegg.v4.0.symbols.gmt [KEGG]
%       - c2.cp.reactome.v4.0.symbols.gmt [Reactome]
%
% For more options, see:
% http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_GSEA_from
%
%

global WorkFolder

% default values (can be specified as input parameters)
p = inputParser;
addParameter(p,'Randseed',ceil(mod(now,1)*1e5), @isscalar);
addParameter(p,'nperm',1e3, @isscalar);
addParameter(p,'Nplot',0, @isscalar);
addParameter(p,'set_min',15, @isscalar);
addParameter(p,'set_max',500, @isscalar);
addParameter(p,'scoring_scheme','weighted', @(x) ismember(x, ...
    {'weighted' 'classic' 'weighted_p2' 'weighted_p1.5'}));
addParameter(p,'Outputfolder', '', @ischar);
addParameter(p,'Outputname', '', @ischar);
addParameter(p,'label', '', @ischar);
addParameter(p,'verbatim', false, @islogical);
addParameter(p,'GSEAfolder', [WorkFolder 'GSEA_java' filesep], @ischar);

parse(p,varargin{:});
p = p.Results;
Outputfolder = p.Outputfolder;
label = p.label;
GSEAfolder = p.GSEAfolder;

% create a temp folder based on the day (will be anyway created by GSEA
tempfolder = ['.' filesep ...
    lower(char(datetime('now','format','MMMdd'))) filesep];
if ~exist(tempfolder,'dir')
    mkdir(tempfolder)
else
    rmdir(tempfolder,'s')
    mkdir(tempfolder)
end
if isempty(Outputfolder)
    Outputfolder = tempfolder;
    if p.Nplot>0
        warnprintf('To plot results, sepcify an ''Outpufolder'' -> Nplot=0')
    end
    p.Nplot = 0;
elseif ~exist(Outputfolder,'dir')
    mkdir(Outputfolder)
end

Inputfile = [tempfolder 'temp.rnk'];
if ischar(Genelist)
    assert(exist(Genelist,'file')>0, 'Input file not found')
    Inputfile = Genelist;
elseif istable(Genelist)
    table2tsv(Genelist, Inputfile, 0);
elseif iscell(Genelist) && size(Genelist,2)<=2
    cell2tsv(Inputfile, Genelist);
else
    error('unknown format for Input')
end


Setfile = [tempfolder 'tempset.gmt'];
if ischar(Geneset)
    switch Geneset
        case 'GOmf'
            Setfile = [GSEAfolder 'c5.mf.v4.0.symbols.gmt'];
        case 'GObp'
            Setfile = [GSEAfolder 'c5.bp.v4.0.symbols.gmt'];
        case 'GOcc'
            Setfile = [GSEAfolder 'c5.cc.v4.0.symbols.gmt'];
        case 'allGO'
            Setfile = [GSEAfolder 'c5.all.v4.0.symbols.gmt'];
        case 'KEGG'
            Setfile = [GSEAfolder 'c2.cp.kegg.v4.0.symbols.gmt'];
        case 'Reactome'
            Setfile = [GSEAfolder 'c2.cp.reactome.v4.0.symbols.gmt'];
        case 'Biocarta'
            Setfile = [GSEAfolder 'c2.cp.biocarta.v4.0.symbols.gmt'];
        case 'CGN'  % cancer gene neighborhoods
            Setfile = [GSEAfolder 'c4.cgn.v4.0.symbols.gmt'];
        case 'CGP'  % chemical and genetic perturbations
            Setfile = [GSEAfolder 'c2.cgp.v4.0.symbols.gmt'];
        case 'BreastCGP'  % chemical and genetic perturbations
            Setfile = [GSEAfolder 'cgp_breast_cancer.gmt'];
        case 'OncoSig'  % oncogenic signatures gene sets
            Setfile = [GSEAfolder 'c6.all.v4.0.symbols.gmt'];
        otherwise
            assert(exist(Geneset,'file')==2, 'Geneset file not found')
            [~,~,ext] = fileparts(Geneset);
            assert(strcmp(ext,'.gmt'),'Geneset file must be a .gmt format')
            Setfile = Geneset;
    end
elseif iscellstr(Geneset)
    assert(size(Geneset,2)==3, 'Geneset as a cellstr should be a nx3 array')
    cell2tsv(Setfile, Geneset);
else
    error('unknown format for Geneset')
end
% check if the file exist locally, otherwise use the
% Broad http address for the file: ftp.broadinstitute.org://pub/gsea/gene_sets/

% add the dataset in the label if not present
if ismember(Geneset, {'GOmf' 'GObp'})
    if isempty(strfind(label, Geneset))
        label = [label '_' Geneset];
    end
    if ~exist(Setfile, 'file')
        [~,Setfile] = fileparts(Setfile);
        Setfile = ['ftp.broadinstitute.org://pub/gsea/gene_sets/' Setfile '.gmt'];
    end
end

%% construct the command
cmd = [ ...
    'java -Xmx1024m -cp ' GSEAfolder 'gsea2-2.1.0.jar' ...
    ' xtools.gsea.GseaPreranked' ...
    ' -gmx ' Setfile ...
    ' -collapse false -mode Max_probe -norm meandiv -nperm ' num2str(p.nperm) ...
    ' -rnk ' Inputfile  ' -scoring_scheme ' p.scoring_scheme ' -rpt_label ' label ...
    ' -include_only_symbols true -make_sets true -plot_top_x ' num2str(p.Nplot) ...
    ' -rnd_seed ' num2str(p.Randseed) ' -set_max ' num2str(p.set_max) ...
    ' -set_min ' num2str(p.set_min) ' -zip_report false' ...
    ' -out ' Outputfolder ' -gui false'];

% and run it
if p.verbatim
    status = system(cmd);
    disp('----------------------------------------------')
    assert(status==0, 'GSEA failed: %s')
else
    [status, output] = system(cmd);
    assert(status==0, 'GSEA failed: %s', output)
end
%% get the right folder

f = dir([Outputfolder filesep label '*GseaPre*']);
f = f([f.isdir]);
Lastrun = f(argmax(cellfun(@datenum,{f.date}))).name;
if p.Nplot>0
    disp(['---> GSEA done Results stored in ' Outputfolder])
    disp([' --->  ' filesep Lastrun])
else
    disp(' ---> GSEA done')
end

%% and get the genes of each set and their scores

% scores are in edb/results.edb
res = XMLElement.parse([Outputfolder filesep Lastrun filesep ...
    'edb' filesep 'results.edb']);

f = dir([Outputfolder filesep Lastrun filesep 'edb' filesep '*rnk']);
Genelist = tsv2cell([Outputfolder filesep Lastrun filesep 'edb' filesep f.name]);
%%
fields = {'GENESET' 'ES' 'NES' 'NP' 'FDR' 'FWER' 'HIT_INDICES' 'ES_PROFILE'};
fieldnames = {'GeneSet' 'Escore' 'NormEscore' 'pval' 'FDR' 'FWER' 'GeneIdx' ...
    'GeneEscore' 'GeneNames' 'LeadEdge'};

values = cell(length(res.children), length(fieldnames));
for i=1:length(res.children)
    for j=1:length(fields)
        values{i,j} = res.children(i).get(fields{j});
    end

    values{i,strcmp(fieldnames, 'GeneNames')} =  ...
        Genelist(1+str2num(values{i,strcmp(fields, {'HIT_INDICES'})}),1)';
end

idx = ismember(fields, {'ES' 'NES' 'NP' 'FDR' 'FWER'});
values(:,idx) = cellfun2(@str2num,values(:,idx));

idx = strcmp(fields, 'GENESET');
values(:,idx) = regexpcellsplit(values(:,idx), '#', 2);

idx = strcmp(fields, 'ES_PROFILE');
values(:,idx) = cellfun2(@str2num,values(:,idx));
idx = strcmp(fields, 'HIT_INDICES');
values(:,idx) = cellfun2(@(x) 1+str2num(x), values(:,idx));

t_results = cell2table(values, 'variablenames', fieldnames);

if isnumeric(t_results.GeneEscore)
    % convert to the cell array of vectors instead of a matrix; this can
    % happen when creating the table if if one row only or all rows have
    % the same number of genes
    t_results.GeneEscore = rowfun_cells(@(x) horzcat(x{:}), ...
        num2cell(t_results.GeneEscore));
    t_results.GeneIdx = rowfun_cells(@(x) horzcat(x{:}), ...
        num2cell(t_results.GeneIdx));
    t_results.GeneNames = rowfun_cells(@(x) horzcat(x{:}), ...
        num2cell(t_results.GeneNames));
end
for i=1:height(t_results)
    if t_results.Escore(i)>0
        t_results.LeadEdge{i} = t_results.GeneNames{i}(1:find(...
            t_results.GeneEscore{i}==max(t_results.GeneEscore{i})));
    else
        t_results.LeadEdge{i} = t_results.GeneNames{i}(find(...
            t_results.GeneEscore{i}==min(t_results.GeneEscore{i})):end);
    end
end

%% cleaning up
if exist(tempfolder,'dir')
    disp('cleaning stuff')
    rmdir(tempfolder,'s')
end
if p.Nplot>0
    if ~isempty(p.Outputname)
        if exist([Outputfolder filesep p.Outputname],'dir')
            if ~rmdir([Outputfolder filesep p.Outputname])
                cnt = 1;
                while exist([Outputfolder filesep p.Outputname '_' num2str(cnt)],'dir') && ...
                        ~rmdir([Outputfolder filesep p.Outputname '_' num2str(cnt)])
                    cnt = cnt+1;
                end
                p.Outputname = [p.Outputname '_' num2str(cnt)];
            end
        end
        movefile([Outputfolder filesep Lastrun], [Outputfolder filesep p.Outputname])
        Lastrun = p.Outputname;
    end
    disp(['Results saved in ' Outputfolder filesep p.Outputname])
    disp('Opening report')
    web([Outputfolder filesep Lastrun filesep 'index.html'],'-browser')
end
