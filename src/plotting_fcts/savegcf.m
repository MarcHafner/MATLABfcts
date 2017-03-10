function savegcf(filename)
% savegcf(filename)

set(gcf,'Renderer','painters')
if nargin==0
    saveas(gcf,get(gcf,'FileName'))
else
    saveas(gcf,filename)
end
