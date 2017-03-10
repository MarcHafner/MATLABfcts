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
    if ismac
        args = '';
    else
        args = 'cat';
    end
end





if exist(outputfile,'file')
    delete(outputfile)
end

if ismac
    inputfiles = strjoin(regexp(inputfiles,' ','split'), '\\ ');
    link = '"/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py"';
    cmd = sprintf('%s %s -o "%s" %s', link, args, outputfile, inputfiles);
    [status, output] = system(cmd);
    assert(status==0, 'combine PDF was not successful:\n\n%s', output)
    disp(output)
else
    d = pwd;
    % absolute path
    inputfiles = [d filesep ReplaceName(inputfiles, '/\', filesep)];
    outputfile = [d filesep ReplaceName(outputfile, '/\', filesep)];
    link = 'pdftk';
    cmd = sprintf('%s "%s" cat %s output "%s"', link, inputfiles, args, outputfile);
    
    [status, output] = system(cmd);
    assert(status==0, 'pdftk was not successful:\n\n%s', output)
end


if exist('del','var') && del
    n = dir(inputfiles);
    [~,f,e] = fileparts(outputfile);
    n = setdiff({n.name}, [f e]);
    n = strcat([fileparts(inputfiles) filesep], n);
    delete(n{:})
end
