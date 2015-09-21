function pathlist = getpath_nogit(folder)

temp = genpath(folder);
temp = regexp(temp,';','split');
temp = temp(cellfun(@isempty,strfind(temp,'.git')));

pathlist = temp{1};
for i=2:length(temp)
    pathlist = [pathlist ';' temp{i}];
end
