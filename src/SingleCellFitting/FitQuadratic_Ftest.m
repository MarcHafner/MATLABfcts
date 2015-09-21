function Fitted_data = FitQuadratic_Ftest(Process_data, plotting)


if ~exist('plotting','var') || isempty(plotting)
    plotting = false;
end
temp_plot='y';

Fitted_data = Process_data;

idx_start = floor(Process_data.parameters.ProcesswindowSize/2);
fit_start = Process_data.T(idx_start);

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0, -30-fit_start],...
    'Upper',[.01, Process_data.parameters.maxT-30-fit_start],...
    'Startpoint',[.002 30]);
% fit with a t^2 of coeeficient  a  from time point  b  with time shift of  c
Fitted_data.fitfct = fittype('(a*(x-b).^2)*(x>b)','options',s);


p_thres = 0.05; % threshold for the F-test against a flat model
windowSize = Process_data.parameters.ProcesswindowSize;

% warning('********************')
for iW=1:Process_data.Nwells

    fprintf('\n** Analyzing %i cells in well %s (%s) **\n', ...
        Process_data.NgoodTraj(iW), Process_data.Wells{iW,3}, ...
        Process_data.WellLabel{iW})


    Fitted_data.fit(iW).r2 = NaN(1,length(Process_data.stats(iW).Endidx));
    Fitted_data.fit(iW).k = Fitted_data.fit(iW).r2;
    Fitted_data.fit(iW).Delay = Fitted_data.fit(iW).r2;
    Fitted_data.fit(iW).endfit = Fitted_data.fit(iW).r2;
    Fitted_data.fit(iW).badidx = Fitted_data.fit(iW).r2;
    Fitted_data.fit(iW).ActivityEndFit = Fitted_data.fit(iW).r2;
    Fitted_data.fit(iW).MaxFitActivity = Fitted_data.fit(iW).r2;
    Fitted_data.fit(iW).flat_fit = false(1,length(Process_data.stats(iW).Endidx));
    Fitted_data.fit(iW).FitToMax = true(1,length(Process_data.stats(iW).Endidx));

    for j=1:Process_data.NgoodTraj(iW)
        if ~Process_data.stats(iW).FRETSurviving(j)
            endfit = Process_data.stats(iW).FRETPreMompTimeIdx(j);
            endfit1 = endfit;
            endfit_max = Process_data.stats(iW).FRETPreMompTimeIdx(j);
        else
            endfit = Process_data.stats(iW).MaxActivityIdx(j);
            endfit1 = endfit;
            endfit_max = min(endfit+2*windowSize, ...
                find(~isnan(Process_data.goodTraj(iW).ICTraj(:,j)),1,'last'));
        end

        traj = Process_data.goodTraj(iW).ICTraj(1:endfit,j);

        if isnan(endfit)
            Fitted_data.fit(iW).badidx(j) = true;
            warning('bad fit at well %s, track %i', Process_data.Wells{iW,3}, j)
            continue
        end




        fprintf('\t % 3i: ', j);
        if endfit1>idx_start

            % redefine the parameter range and seed
            [MinIC, MinTimeIC] = nanmin(traj);
            MinTimeIC = Process_data.T(MinTimeIC);
            startpt = [ (traj(endfit)-MinIC)/ (Process_data.T(endfit)-MinTimeIC)^2 ...
                MinTimeIC];

            set(s,'Startpoint',startpt,'upper', ...
                [.01, Process_data.T(endfit)-30-fit_start] )
            f = fittype('(a*(x-b).^2)*(x>b)','options',s);
            [coeff, gof] = fit(Process_data.T(idx_start:endfit)'-fit_start, ...
                traj(idx_start:endfit), f);

            r2 = gof.rsquare;
            p = Ftest_flat(traj(idx_start:endfit), gof);
            doing_peaks = (r2<.5) || (p>p_thres);
            if doing_peaks
                fprintf('fit(r=%5.2f, p=%.2f) -> peaks: ', gof.rsquare, p)
            end
        else
            fprintf('too early max -> peaks: ');
            endfit = nanmin([endfit_max, ...
                find(Process_data.goodTraj(iW).ICTraj(:,j)==...
                max(Process_data.goodTraj(iW).ICTraj(:,j)))+2*windowSize]);
            doing_peaks = true;
            r2 = 0;
        end


        if doing_peaks
            traj = smooth(Process_data.goodTraj(iW).C8Activity(...
                min(windowSize,endfit_max-2):min(endfit_max, ...
                length(Process_data.T)-2*windowSize),j),3);
            if mean(diff(traj)>0)>.9
                idx2 = argmax(traj);
                fprintf('special case of a short increase trajectory')
            else
                [~, idx2] = findpeaks(traj,'MINPEAKHEIGHT', Process_data.stats(iW).MaxActivity(j)/3);
            end

            idx2 = idx2+windowSize-1;
            idx2 = idx2( idx2>idx_start);
            %                 if isempty(idx2) && any(smooth(...
            %                         Process_data.goodTraj(iW).C8Activity(1:min(endfit_max, ...
            %                         length(Process_data.T)-2*windowSize),j),3) > 1e-3)
            %                     [~, idx2] = findpeaks(smooth(...
            %                         Process_data.goodTraj(iW).C8Activity(1:min(endfit_max, ...
            %                         length(Process_data.T)-2*windowSize),j),3) , ...
            %                         'SORTSTR', 'descend','MINPEAKHEIGHT', 1e-3);
            %                     idx2 = idx2( idx2>idx_start);
            %                 end
            if isempty(idx2)
                Fitted_data.fit(iW).flat_fit(j) = true;
                fprintf(' flat fit: r=%4.2f, k=0\n', r2);
            else

                r2s = 0*idx2;
                ps = r2s;
                c2s = {};

                for i2=1:length(idx2)
                    endfit2 = min(idx2(i2), ...
                        find(~isnan(Process_data.goodTraj(iW).ICTraj(:,j)),1,'last'));
                    traj = Process_data.goodTraj(iW).ICTraj(1:endfit2,j);

                    [MinIC, MinTimeIC] = nanmin(traj);
                    MinTimeIC = Process_data.T(MinTimeIC);
                    startpt = [ (traj(endfit2)-MinIC)/ (Process_data.T(endfit2)-MinTimeIC)^2 ...
                        MinTimeIC];

                    set(s,'Startpoint',startpt,'upper',[.01, Process_data.T(endfit2)-30-fit_start] )
                    f = fittype('(a*(x-b).^2)*(x>b)','options',s);
                    [c2s{i2}, gof2] = fit(Process_data.T(idx_start:endfit2)'-fit_start, ...
                        traj(idx_start:endfit2),f);
                    r2s(i2) = gof2.rsquare;

                    % perform the F-test against flat line:
                    ps(i2) = Ftest_flat(traj(idx_start:endfit2), gof2);

                end
                disp([r2s ps])
                [r2, i2] = max(r2s.*(ps<.05));
                coeff2 = c2s{i2};
                endfit2 = idx2(i2);

                if all( ps>p_thres )
                    Fitted_data.fit(iW).flat_fit(j) = true;
                    fprintf(' flat fit: r=%4.2f, k=0\n', r2);

                else
                    Fitted_data.fit(iW).FitToMax(j) = endfit == endfit2;
                    Fitted_data.fit(iW).r2(j) = r2;
                    Fitted_data.fit(iW).k(j) = coeff2.a;
                    Fitted_data.fit(iW).Delay(j) = coeff2.b+fit_start;
                    Fitted_data.fit(iW).endfit(j) = Process_data.T(endfit2);
                    Fitted_data.fit(iW).ActivityEndFit(j) = ...
                        Process_data.goodTraj(iW).C8Activity(endfit2,j);
                    fprintf('; best fit at T=%3i: r=%4.2f, k=%6.3e, p=%.3f\n',....
                        Process_data.T(endfit2), r2, coeff2.a, ps(i2));
                end
            end
        else

            Fitted_data.fit(iW).r2(j) = gof.rsquare;
            Fitted_data.fit(iW).k(j) = coeff.a;
            Fitted_data.fit(iW).Delay(j) = coeff.b+fit_start;
            Fitted_data.fit(iW).endfit(j) = Process_data.T(endfit);
            Fitted_data.fit(iW).ActivityEndFit(j) = ...
                Process_data.goodTraj(iW).C8Activity(endfit,j);
            fprintf(' --> good fit: r=%4.2f, k=%6.3e, tau=%.0f, p=%.3f <--\n', ...
                gof.rsquare, coeff.a, coeff.b+fit_start, p)
        end

        if Fitted_data.fit(iW).flat_fit(j)
            Fitted_data.fit(iW).MaxFitActivity(j) = 0;
            Fitted_data.fit(iW).r2(j) = 0;
            Fitted_data.fit(iW).k(j) = 0;
            Fitted_data.fit(iW).Delay(j) = 0;
            Fitted_data.fit(iW).endfit(j) = Process_data.T(endfit);
            Fitted_data.fit(iW).ActivityEndFit(j) = 0;
        else
            Fitted_data.fit(iW).MaxFitActivity(j) = ...
                2*Fitted_data.fit(iW).k(j)*(Fitted_data.fit(iW).endfit(j)-...
                Fitted_data.fit(iW).Delay(j));
        end


        if (abs(plotting)>=1 && j<abs(plotting) && temp_plot=='y') || ...
                (abs(plotting)<1 && plotting~=0 && Fitted_data.fit(iW).r2(j)<abs(plotting) && temp_plot=='y') || ...
                (plotting~=0 && doing_peaks && temp_plot=='y')
            figure(995);clf
            subplot(1,3,[1 2])
            hold on
            idx = find(~isnan(Process_data.goodTraj(iW).ICTraj(:,j)),1,'last');
            plot(Process_data.T, Process_data.goodTraj(iW).FRET(:,j),'c','linewidth',1)
            plot(Process_data.T(1:idx), Process_data.goodTraj(iW).ICTraj(1:idx,j),'m','linewidth',2)
            plot(Process_data.stats(iW).MaxActivityTime(j)*[1 1], ...
                [0 Process_data.goodTraj(iW).ICTraj(Process_data.stats(iW).MaxActivityIdx(j),j)], 'ok-','linewidth',1)
            if Process_data.stats(iW).FRETPreMompTime(j)>0
                plot(Process_data.stats(iW).FRETPreMompTime(j)*[1 1], ...
                    [0 Process_data.stats(iW).MompThres(j)], '.g-')
            end
            if endfit1>idx_start
                plot(Process_data.T(idx_start:endfit1), coeff(Process_data.T(idx_start:endfit1)-fit_start),'b','linewidth',2)
            end
            if doing_peaks && ~Fitted_data.fit(iW).flat_fit(j)
                plot(Process_data.T(idx_start:endfit2), coeff2( Process_data.T(idx_start:endfit2)' -fit_start),'r','linewidth',2)
                plot(Process_data.T(endfit2)*[1 1], [0 Process_data.goodTraj(iW).ICTraj(endfit2,j)], 'ok-','linewidth',1)
            elseif Fitted_data.fit(iW).flat_fit(j)
                 plot(Process_data.T([1 idx]), [0 0],'r','linewidth',2)
            end
            plot(Fitted_data.fit(iW).endfit(j)+[-30 0 10], ...
                Fitted_data.fit(iW).k(j)*(Fitted_data.fit(iW).endfit(j)-...
                Fitted_data.fit(iW).Delay(j))^2 + ...
                [-30 0 10]*Fitted_data.fit(iW).MaxFitActivity(j), 'xk-','markersize',15,'linewidth',2)
            xlim(Process_data.T([1 end]))

            subplot(133)
            hold on
            plot(Process_data.T(1:idx), [0;diff(Process_data.goodTraj(iW).ICTraj(1:idx,j))],'k','linewidth',1)
            plot(Process_data.T(1:idx), Process_data.goodTraj(iW).C8Activity(1:idx,j),'b','linewidth',2)
            plot(Process_data.T([1 idx]), [0 0],'c')
            plot(Process_data.stats(iW).MaxActivityTime(j)*[1 1], ...
                [0 Process_data.stats(iW).MaxActivity(j)], 'ok-')
            if doing_peaks  && ~Fitted_data.fit(iW).flat_fit(j)
                plot(Process_data.T(endfit2)*[1 1], [0 Process_data.goodTraj(iW).C8Activity(endfit2,j)], 'ok-','linewidth',1)
            end
            if ~isnan(Process_data.stats(iW).FRETPreMompTimeIdx(j))
                plot(Process_data.stats(iW).FRETPreMompTime(j)*[1 1], ...
                    [0 Process_data.stats(iW).FRETMompActivity(j)], '-r.')
            end
            plot(Fitted_data.fit(iW).endfit(j), Fitted_data.fit(iW).MaxFitActivity(j), 'xk','markersize',15)
            xlim(Process_data.T([1 end]))
            if plotting>0
                temp_plot = input('    Do you want more? y/n [y]: ', 's');
                if isempty(temp_plot)
                    temp_plot = 'y';
                end
            else
                pause(.5)
            end
        elseif abs(plotting)>0 && temp_plot~='y'
            break
        else
            continue
        end
    end
    temp_plot = 'y';
end
end


function p = Ftest_flat(traj, gof)
RSS_quad = gof.sse;
RSS_flat = sum((traj-mean(traj)).^2);
F = ((RSS_flat-RSS_quad)/1)/(RSS_quad/(length(traj)-2));
p = 1-fcdf(F, 1, length(traj)-2);
end
