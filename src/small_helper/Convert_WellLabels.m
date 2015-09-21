function welllables = Convert_WellLabels(labels)

welllables = cellfun2(@(x) [x(1) num2str(str2double(x(2:end)),'%02i')],labels);
