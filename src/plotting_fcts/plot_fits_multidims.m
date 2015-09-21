function a = plot_fits_multidims(t_fits, varargin)

%  a = plot_fits_multidims(t_fits, varargin)
%       inputs:
%       - xplotkey
%       - yplotkey
%       - xaxisval      (concentration; default is range of data)
%       - ydata         (IC or GI)
%       - colorkey
%       - raw_data
%       - xtransform
%       - ytransform
%       - axischanges
%       - xspacing
%       - yspacing
%       - yval_lines
%

Generate_Plotting_parameters

noyplot = 'noyplot';
nocplot = 'nocplot';

p = inputParser;
addParameter(p, 'xplotkey', @isstr)
addParameter(p, 'yplotkey', noyplot, @isstr)
addParameter(p, 'xaxisval', [], @isnumeric)
addParameter(p, 'ydata', [], @isstr)
addParameter(p, 'colorkey', nocplot, @isstr)
addParameter(p, 'raw_data', true, @boolean)
addParameter(p, 'xtransform', @log10, @(x)isa(x,'function_handle'));
addParameter(p, 'ytransform', @(x)x, @(x)isa(x,'function_handle'))
addParameter(p, 'axischanges', @(x) set(x,'fontsize',6), @(x) isa(x,'function_handle'))
addParameter(p, 'plotcolors', Plotting_parameters.colors, @isnumeric)
addParameter(p, 'xspacing', .03, @isnumeric)
addParameter(p, 'yspacing', .07, @isnumeric)
addParameter(p, 'yval_lines', [0 1], @isnumeric)

parse(p,varargin{:})
p = p.Results;

xplotkeys = unique(t_fits.(p.xplotkey));
if strcmp(p.yplotkey, noyplot)
    yplotkeys = 1;
    t_fits = [t_fits table(ones(height(t_fits),1), 'variablenames', {noyplot})];
else
    yplotkeys = unique(t_fits.(p.yplotkey));
end
if strcmp(p.colorkey, nocplot)
    colorkeys = 1;
    t_fits = [t_fits table(ones(height(t_fits),1), 'variablenames', {nocplot})];
else
    colorkeys = unique(t_fits.(p.colorkey));
end
if strcmpi(p.ydata,'IC') || ~isvariable(t_fits, 'nGI_fit')
    ydata = 'Relative count';
elseif strcmpi(p.ydata,'GI')
    ydata = 'Normalized relative growth';
else
    ydata = 'Normalized growth inhibition';
end

%
t_fits = TableToCategorical(t_fits, {p.xplotkey, p.yplotkey, p.colorkey});

%%

nCols = length(xplotkeys);
nRows = length(yplotkeys);

xspacing = p.xspacing;
axis_width = (1-(nCols+1)*xspacing)/nCols;
yspacing = p.yspacing;
axis_height = (1-(nRows+1)*yspacing)/nRows;

for ixp=1:nCols
    for iyp = 1:nRows

        a(iyp, ixp) = get_newaxes([xspacing*1.5+(ixp-1)*(axis_width+xspacing) ...
            yspacing*1.5+(nRows-iyp)*(axis_height+yspacing) axis_width axis_height],1,...
            'fontsize',6);

        title([p.xplotkey '=' AnyToString(xplotkeys(ixp)) '; ' ...
            p.yplotkey '=' AnyToString(yplotkeys(iyp))])

        if isempty(p.xaxisval)
            idx = t_fits.(p.xplotkey)==xplotkeys(ixp) & ...
                t_fits.(p.yplotkey)==yplotkeys(iyp);
            if ~any(idx), continue, end
            xvals = t_fits.Conc{idx};
        else
            xvals = p.xaxisval;
        end
        xvals = p.xtransform(xvals);
        for i=1:length(p.yval_lines)
            plot([min(xvals(:)) max(xvals(:))], p.yval_lines(i)*[1 1], ...
                '-', 'color', [.8 .8 .8])
        end

        for iC = 1:length(colorkeys)
            subt = t_fits(t_fits.(p.xplotkey)==xplotkeys(ixp) & ...
                t_fits.(p.yplotkey)==yplotkeys(iyp) & ...
                t_fits.(p.colorkey)==colorkeys(iC),:);

            for i=1:height(subt)
                if isempty(p.xaxisval)
                    xvals = 10.^(log10(min(subt.Conc{i})):.01:log10(max(subt.Conc{i})));
                else
                    xvals = p.xaxisval;
                end

                Conc = p.xtransform(subt.Conc{i});
                if strcmp(ydata, 'Relative count')
                    h(iC) = plot(p.xtransform(xvals), p.ytransform(subt.fit{i}(xvals)), '-', ...
                        'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    if p.raw_data
                        plot(Conc, subt.RelCellCnt{i}, '.', ...
                            'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    end
                elseif strcmp(ydata, 'Normalized relative growth')
                    h(iC) = plot(p.xtransform(xvals), p.ytransform(subt.nGI_fit{i}(xvals)), '-', ...
                        'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    if p.raw_data
                        plot(Conc, p.ytransform(subt.nRelGrowth{i}), '.', ...
                            'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    end
                else
                    h(iC) = plot(p.xtransform(xvals), p.ytransform(subt.GI_fit{i}(xvals)), '-', ...
                        'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    if p.raw_data
                        plot(Conc, p.ytransform(subt.RelGrowth{i}), '.', ...
                            'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    end
                end
            end
        end

        p.axischanges(gca)
        if iyp==nRows, xlabel(gca,'Concentration (\muM)','fontweight','bold'), end
        if ixp==1, ylabel(gca, ydata,'fontweight','bold'), end

        if ixp==nCols && iyp==nRows
            hl = legend(h, strcat(p.colorkey, '=', AnyToString(colorkeys)), ...
                'fontsize',6, 'orientation', 'horizontal');
            set(hl, 'position', [.01 .005 .98 .03])
        end
    end
end

if nargout==0
    clear a
end
