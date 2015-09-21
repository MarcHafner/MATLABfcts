function [MSE, PredOutputs, nVar, B, q2, MCC, details] = LV1outCV_ElasticNet(Inputs, Outputs, opt)
% [MSE, PredOutputs, nVar, B, q2, MCC, details] = ...
%       LV1outCV_ElasticNet(Inputs, Outputs, opt)
%
%       opt = {'lambdas' 'alphas'}


nCL = length(Outputs);
lambdas = 10.^[-4:.15:3];
alphas = 3.^[-4:.2:0];


if exist('opt','var')
    vars = {'lambdas' 'alphas'};
    for i=1:length(vars)
        if isfield(opt,vars{i})
            disp(vars{i});
            eval([vars{i}  ' = opt.' vars{i} ';'])
        end
    end
end


allMSE = 0*alphas'*lambdas;
allMCC = allMSE;
nVars = allMSE;

alltempB = zeros(nCL, size(Inputs,2), length(lambdas), length(alphas));
allPredOutputs = zeros(nCL,length(lambdas),length(alphas));

[BinOutputs, BinarizationThres] = UnsatMedianBinarizeColumn_v2(Outputs);
tic
fprintf([num2str(nCL) ' Cell lines: ']);

nCore = intmax;
if length(alphas)<5
    nCore = 0;
end
parfor (i = 1:length(alphas), nCore )
% for i = 1:length(alphas)

    tempnVar = zeros(nCL,length(lambdas));
    Prediction = zeros(nCL,length(lambdas));
    tB = zeros(nCL, size(Inputs,2), length(lambdas));

    for j=1:nCL
        testidx = setdiff(1:nCL,j);

        shiftInputs = mean(Inputs(testidx,:));
        shiftOutputs = mean(Outputs(testidx));

        NormInput = Inputs(testidx,:)-repmat(shiftInputs,nCL-1,1);
        NormOutputs = Outputs(testidx)-shiftOutputs;

        [B, FitInfo] = lasso(NormInput,NormOutputs,...
            'alpha',alphas(i),'Standardize', false,'Lambda',lambdas);

        Prediction(j,:) = (Inputs(j,:)-shiftInputs)*B + ...
            FitInfo.Intercept(:)' + shiftOutputs;

         tB(j,:,:) = B;

        tempnVar(j,:) = FitInfo.DF;
    end

    alltempB(:,:,:,i) = tB;
    allPredOutputs(:,:,i) = Prediction;
    allMSE(i,:) = sum((Prediction - repmat(Outputs,1,length(lambdas))).^2);

    tempMCC = 0*lambdas;
    for j=1:length(lambdas)
        tempMCC(j) = PredictionScoring(Prediction(:,j)>BinarizationThres,BinOutputs);
    end
    allMCC(i,:) = tempMCC;

    nVars(i,:) = mean(tempnVar);
    fprintf([' ' num2str(i)]);

end
fprintf('\n');
toc
%%

mSE = min(allMSE(:));

MCC = max(allMCC(:) - 2*(allMSE(:)>(mSE*1.02)));
[alpha, lambda] = find( (nVars.*(allMCC==MCC)) == min(min(nVars(allMCC(:)==MCC))),1,'first');

nVar = nVars(alpha,lambda);
MSE = allMSE(alpha,lambda);
PredOutputs = allPredOutputs(:,lambda,alpha);
B = (alltempB(:,:,lambda,alpha));

SSquares = sum( (Outputs-mean(Outputs)).^2 );
q2 = 1 - (MSE/SSquares);

if nargout == 7
    details.lambdas = lambdas;
    details.alphas = alphas;
    details.allMSE = allMSE;
    details.allPredOutputs = allPredOutputs;
    details.nVars = nVars;
    details.alltempB = alltempB;

    details.Falpha = alphas(alpha);
    details.Flambda = lambdas(lambda);

    NormInput = Inputs-repmat(mean(Inputs),nCL,1);
    NormOutputs = Outputs-mean(Outputs);

    [details.FB, FitInfo] = lasso(NormInput,NormOutputs,...
        'alpha',details.Falpha,'Standardize',false,'Lambda',details.Flambda);

    details.FPrediction = NormInput*details.FB + FitInfo.Intercept' + mean(Outputs);

end
