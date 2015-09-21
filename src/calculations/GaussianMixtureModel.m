function [bestfit, bestIC] = GaussianMixtureModel(data,varargin)
% function [bestfit, bestIC] = GaussianMixtureModel(data,varargin)
%
%       varargin:
%           - plotting  [false]
%           - useAIC    [false]
%           - nGaussMax [6]
%           - RegularizeValue [.0002*range(data)]
%
%



%% inputs and parser
p = inputParser;

addParameter(p,'plotting',false, @islogical);
addParameter(p,'useAIC',false,  @islogical);
addParameter(p,'nGaussMax',6, @isnumeric);
addParameter(p,'RegularizeValue',[], @isnumeric);

parse(p, varargin{:})
p = p.Results;
for i = {'plotting' 'useAIC' 'nGaussMax' 'RegularizeValue'}
    eval([i{:} ' = p.' i{:} ';'])
end




delta = max(data)-min(data);

if isempty(RegularizeValue)
    RegularizeValue = delta*.0002;
end

if isrow(data);data=data';end


usedIC = zeros(nGaussMax,1);
bestIC = Inf;
Ngauss = 0;

warning('off','stats:gmdistribution:FailedToConverge')
warning('off','stats:gmdistribution:FailedToConvergeReps')

for i=1:nGaussMax
    obj = gmdistribution.fit(data,i,'replicates',50,'Regularize',RegularizeValue);
    if useAIC
        usedIC(i) = obj.AIC;
    else
        usedIC(i) = obj.BIC;
    end

    if plotting
        disp(['Iteration ' num2str(i) ': AIC=' num2str(obj.AIC) ' ; BIC=' num2str(obj.BIC)])
    end

    if i>0 && usedIC(i)<bestIC
        bestIC = usedIC(i);
        bestfit = obj;
        Ngauss = i;
    end

    if i>3 && all(usedIC(i)>usedIC(i-[1:2]))
        break
    end
end
%%

warning('on','stats:gmdistribution:FailedToConverge')
warning('on','stats:gmdistribution:FailedToConvergeReps')

if plotting
    washold = ishold;


    X1 = (min(data)-.3*delta):(delta/max(10,sqrt(length(data)))):(max(data)+.3*delta);
    X = X1;
    for i=1:Ngauss
        X = [X bestfit.mu(i)+[(-.5*delta):(.2*delta/sqrt(length(data))):(.5*delta)]*bestfit.Sigma(:,:,i)^.5];
    end
    X = sort(X);
    %%

    n = hist(data,X1);
    Y = pdf(bestfit,X');
    h(1) = plot(X,max(n)*Y/max(Y)/sum(n),'-b','linewidth',3);
    hold on

    for i=1:length(bestfit.mu)
        Y1 = pdf('norm',X,bestfit.mu(i),bestfit.Sigma(:,:,i)^.5)*bestfit.PComponents(i);
        h(2+i) = plot(X,max(n)*Y1/max(Y)/sum(n),'r--','linewidth',2);
    end

    h(2) = bar(X1,n/sum(n),'facecolor','none','linewidth',2);

    if ~washold
        hold off
    end

    legend(h(1:3),{'Gaussian mixture model', 'Cell line distribution', 'individual models'},...
        'location','northwest','fontweight','bold','fontsize',8)
    xlabel('data','fontweight','bold','fontsize',10)
    ylabel('Probability distribution','fontweight','bold','fontsize',10)
    set(gca,'fontweight','bold','fontsize',8)

end
