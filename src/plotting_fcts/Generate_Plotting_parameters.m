

global Plotting_parameters
if isempty(Plotting_parameters)

    Plotting_parameters.cmapRB = [ [.9:.001:1]' [0:.01:1]' [0:.01:1]' ;
        [1:-.01:0]' [1:-.01:0]' [1:-.001:.9]' ].^1.2;
    Plotting_parameters.cmapBR = Plotting_parameters.cmapRB(end:-1:1,:);

    Plotting_parameters.cmapRW = [ [.9:.001:1]' [0:.01:1]' [0:.01:1]' ];
    Plotting_parameters.cmapWR = Plotting_parameters.cmapRW(end:-1:1,:);

    Plotting_parameters.cmapBW = [ [0:.01:1]' [0:.01:1]' [.9:.001:1]' ];
    Plotting_parameters.cmapWB = Plotting_parameters.cmapBW(end:-1:1,:);

    Plotting_parameters.cmapYBB = [ [1:-.009:.1]' [1:-.009:.1]' [.01:-.0001:0]' ;
        [0:.001:.1]' [0:.001:.1]' [0:.01:1]' ].^.7;
    Plotting_parameters.cmapBBY = Plotting_parameters.cmapYBB(end:-1:1,:);

    Plotting_parameters.cmapYGyBl = .2 +[
        [.8:-.005:.3]' [.8:-.005:.3]' [0:.003:.3]'
        [.3:-.002:.1]' [.3:-.002:.1]' [.3:.005:.8]'
        ].^1.5;

    Plotting_parameters.cmapBkY = [ [0:.006:.6]' [0:.006:.6]' [0:.0005:.05]' ;
        [.6:.004:1]' [.6:.004:1]' [.05:.0015:.2]' ].^1.2;


    Plotting_parameters.cmapGrWP = [ [0:.01:1]' [.3:.007:1]' [.2:.008:1]' ;
        [1:-.005:.5]' [1:-.0095:.05]' [1:-.004:.6]' ].^1.1;

    Plotting_parameters.cmapWP = (1.2*(Plotting_parameters.cmapGrWP(end/2:end,:)-.7))+.75;
    Plotting_parameters.cmapWP(Plotting_parameters.cmapWP>1) = 1;
    Plotting_parameters.cmapWP(Plotting_parameters.cmapWP<0) = 0;


    Plotting_parameters.cmapGrWBr = [ [0:.01:1]' [.6:.004:1]' [.4:.006:1]'
        [1:-.0035:.65]' [1:-.0055:.45]' [1:-.0095:.05]' ].^1.1;
    Plotting_parameters.cmapBrWGr = Plotting_parameters.cmapGrWBr(end:-1:1,:);

    Plotting_parameters.cmapWBr = Plotting_parameters.cmapGrWBr(end/2:end,:).^1.1;
    % for i=1:3
    %     Plotting_parameters.cmapWBr(:,1) = max([zeros(length(Plotting_parameters.cmapWBr),1) ...
    %         min([ones(length(Plotting_parameters.cmapWBr),1) ...
    %         (1.2*(Plotting_parameters.cmapWBr(:,i)-.75))+.75],[],2) ],[],2);
    Plotting_parameters.cmapWGr = Plotting_parameters.cmapBrWGr(end/2:end,:).^1.1;

    Plotting_parameters.colorsNames = 'brmcgykwbrmcgykw';
    Plotting_parameters.linetypes = {'-' '--' '.-' ':'};
    Plotting_parameters.markers = 'o^sdv>xp<';
    % Plotting_parameters.markersizes=[10.5 10 10 10 10 10 10]*1.1;
    % Plotting_parameters.markersizesScatter=[14 12 13 12 12 12 12 12]*12;
    Plotting_parameters.markersizes=[10.5 10 10 10 10 10 10]/1.1;
    Plotting_parameters.markersizesScatter=[14 12 13 12 12 12 12 12]*8;
    Plotting_parameters.gray = [.7 .7 .7];
    Plotting_parameters.SelectedGrays = [.2;.55;.8]*[1 1 1];

    Plotting_parameters.colors = [
        250 50 50
        35 60 255
        170 95 194
        95 150 100
        249 226 20
        256 256 190
        0 0 0]/256;

    for iColorsInPlotting=1:size(Plotting_parameters.colors,1)
        Plotting_parameters.Lightcolors(iColorsInPlotting,:) = ...
            min([1 1 1;(Plotting_parameters.colors(iColorsInPlotting,:)+.15).^.8]);
        Plotting_parameters.Darkcolors(iColorsInPlotting,:) = ...
            max([0 0 0;(Plotting_parameters.colors(iColorsInPlotting,:)-.15)]).^1.2;
    end
    Plotting_parameters.Lightcolors(3,:) = ...
        min([1 1 1;(Plotting_parameters.colors(3,:)+.05).^.9]);

    Plotting_parameters.ProfilesColorsIdx = [5 4 3 3];

    Plotting_parameters.axes = {'fontsize',6};
    Plotting_parameters.legend = {'fontsize',6,'fontweight','bold'};
    Plotting_parameters.axislabel = {'fontsize',8,'fontweight','bold'};
    Plotting_parameters.title = {'fontsize',10,'fontweight','bold'};
    Plotting_parameters.plotText = {'fontsize' 8 'fontweight' 'bold'};

    Plotting_parameters.logticks = @(x,s,y) 10.^(floor(x):s:ceil(y));
    Plotting_parameters.logticklab = @(x,s,y) num2cellstr(10.^(floor(x):s:ceil(y)),'%5.3g');

    clear iColorsInPlotting
end
