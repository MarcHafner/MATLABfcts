function merge_pdf(inputfiles, outputfile, del, args)
%  merge_pdf(inputfiles, outputfile, del, args)
%
%   relies on the 'server' version of pdftk (see
%   https://www.pdflabs.com/tools/pdftk-server/) pdftk should be accesible
%   from the current folder
%
%   inputfiles  string readable in DOS formet (like ' any*.pdf ')
%   outpfile    name of a pdf file for the output
%   del         delete the original files (leave only the merged file)
%   args        optional arguments for the pdftk command (default is ' cat ')
%
%   WARNING: overwritting the output file without notice
%
%
%

if ~exist('args','var')
    args = 'cat';
end




d = pwd;
% absolute path needed for MS-DOS
inputfiles = [d filesep ReplaceName(inputfiles, '/\', filesep)];
outputfile = [d filesep ReplaceName(outputfile, '/\', filesep)];

if exist(outputfile,'file')
    delete(outputfile)
end

cmd = sprintf('pdftk %s %s output %s', inputfiles, args, outputfile);

[status, output] = system(cmd);
assert(status==0, 'pdftk was not successful:\n\n%s', output)

if exist('del','var') && del
    n = dir(inputfiles);
    [~,f,e] = fileparts(outputfile);
    n = setdiff({n.name}, [f e]);
    n = strcat([fileparts(inputfiles) filesep], n);
    delete(n{:})
end
