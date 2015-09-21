function [q2, CV_Pred, B, r2, fit_output, details] = LV1outCV_nestedPLSRvip(Inputs, Outputs, useVIP)

% function [q2, CV_Pred, B, r2, fit_output, details] = ...
%     LV1outCV_nestedPLSRvip(Inputs, Outputs, useVIP)
%
%   useVIP = 0/1
%%

if ~exist('useVIP','var')
    useVIP = 1;
else
    disp(['using VIP : ' num2str(useVIP)])
end

nCL = size(Outputs,1);
nVar = size(Inputs,2);

if size(Inputs,1) ~= nCL
    error('Wrong number of samples in the Inputs')
end
if size(Outputs,2) ~= 1
    error('Wrong number of columns in the output')
end

individual_r2 = NaN(nCL,1);
individual_q2 = NaN(nCL,1);
CV_Pred = NaN(nCL,1);
individual_models = NaN(nVar,nCL);
steptime = zeros(nCL,1);
vip_cutoff = zeros(nCL,1);

opt.nCore = 0;

totime = tic;

fprintf('%i cell lines', nCL);

parfor i=1:nCL
%     opt = [];for i=1:nCL
    time1 = tic;

    testidx = setdiff(1:nCL,i);

    shiftInputs = mean(Inputs(testidx,:));
    shiftOutputs = mean(Outputs);

    NormInput = Inputs(testidx,:)-repmat(shiftInputs,nCL-1,1);
    NormOutputs = Outputs(testidx)-shiftOutputs;

        if useVIP
            [~, ~, ~, ~, ~, ...
                ~, individual_models(:,i), individual_r2(i), individual_q2(i), details_vip] = ...
                LV1outCV_PLSRvip(NormInput, NormOutputs, opt);
            vip_cutoff(i) = details_vip.VIPcutoff;
        else
            [~, individual_models(:,i), individual_r2(i), individual_q2(i)] = ...
                LV1outCV_PLSRvip(NormInput, NormOutputs, opt);
            vip_cutoff(i) = 0;
        end

    CV_Pred(i) = (Inputs(i,:) - mean(Inputs(testidx,:)))*individual_models(:,i) + mean(Outputs(testidx));

    steptime(i) = toc(time1);
end

NormInput = Inputs-repmat(mean(Inputs),nCL,1);

q2 = 1 - (sum((CV_Pred - Outputs).^2)/sum( (Outputs-mean(Outputs)).^2 ));

fprintf(' ; done in %.2fs, %.1f (av, max); av # vars = %.1f (VIP); q2=%.3f ; ', mean(steptime), toc(totime), mean(sum(individual_models~=0)), q2)

details = [];
if useVIP
        [~, ~, ~, ~, ~, ...
            fit_output, B, r2, ~, temp_details] = ...
            LV1outCV_PLSRvip(NormInput, Outputs);
        details.VIPcutoff = temp_details.VIPcutoff;
    else
        [fit_output, B, r2, ~, temp_details] = ...
            LV1outCV_PLSRvip(NormInput, Outputs);
        details.VIPcutoff = 0;
end


if useVIP
    details.optnComp = temp_details.optnComp;
    details.optnComp2 = temp_details.optnComp2;
    details.VIPcutoff = temp_details.VIPcutoff;
    details.Allvip_cutoff = vip_cutoff;
else
    details.optnComp = temp_details.optnComp;
    details.optnComp2 = 0;
    details.VIPcutoff = 0;
    details.Allvip_cutoff = 0;
end

details.individual_models = individual_models;
details.individual_r2 = individual_r2;
details.individual_q2 = individual_q2;

fprintf('r^2=%.3f . over\n', r2);
