function save_cf(filename)

drawnow

if nargin==1
    saveas(gcf,filename)
else
    saveas(gcf,get(gcf,'FileName'))
end
