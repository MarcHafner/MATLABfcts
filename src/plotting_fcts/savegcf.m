function savegcf(filename)
% savegcf(filename)

if nargin==0
    saveas(gcf,get(gcf,'FileName'))
else
    saveas(gcf,filename)
end
