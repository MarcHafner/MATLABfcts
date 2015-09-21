function [fit_output, B, r2, chosen_q2, details, ...
    fit_output_vip, B_vip, r2_vip, chosen_q2_vip, details_vip] = ...
    LV1outCV_PLSRvip(Inputs, Outputs, opt)
% function [fit_output, B, r2, chosen_q2, details, ...
%     fit_output_vip, B_vip, r2_vip, chosen_q2_vip, details_vip] = ...
%     LV1outCV_PLSRvip(Inputs, Outputs, opt)
%
%   opt :
%         nComps
%         nCore
%         nComps2
%         allVIPthres
%
%
%       doVIP = nargout>4;
%       useVIP = nargout>5;
%
%
%   Marc Hafner
%
%%
doVIP = nargout>4;
useVIP = nargout>5;

nComp = 15;
if abs(nComp)>=rank(Inputs)-2
    nComp = sign(nComp)*max(1,(rank(Inputs)-3));
end

if all(Outputs==mean(Outputs))
    fit_output = Outputs;
    B = zeros(size(Inputs,2),1);
    r2 = 1;
    chosen_q2 = 0;
    details.optnComp = 0;
    details.VIP = B;
    fit_output_vip = Outputs;
    B_vip = zeros(size(Inputs,2),1);
    r2_vip = 1;
    chosen_q2_vip = 0;
    details_vip.optnComp = 0;
    details_vip.VIPcutoff = 0;
    details_vip.optnComp2 = 0;
    return
end

nComps = 1:nComp;
nComps2 = nComps;
allVIPthres = 0:.05:3.5;

nCore = intmax;


if exist('opt','var')
    if isfield(opt,'nComps')
        nComps = opt.nComps;
    end
    if isfield(opt,'nCore')
        nCore = opt.nCore;
    end
    if isfield(opt,'nComps2')
        nComps2 = opt.nComps2;
    end
    if isfield(opt,'allVIPthres')
        allVIPthres = opt.allVIPthres;
    end
end

nCL = size(Outputs,1);
nVar = size(Inputs,2);
nVIP = length(allVIPthres);

allPredOutputs = NaN(nCL,length(nComps),'single');

if useVIP
    allPredOutputs_vip = NaN(nCL,length(nComps),nVIP,length(nComps2),'single');
    nVar_vip = zeros(nCL,length(nComps),nVIP,length(nComps2),'single');
end


if length(nComps)<5
    nCore = 0;
end


parfor (iCL=1:nCL, nCore)
%     for iCL=1:nCL

    testidx = setdiff(1:nCL,iCL);

    shiftInputs = mean(Inputs(testidx,:));
    NormInput = Inputs(testidx,:)-repmat(shiftInputs,nCL-1,1);

    shiftOutputs = mean(Outputs(testidx));
    NormOutputs = Outputs(testidx)-shiftOutputs;

    if all(abs(NormOutputs)<1e-10)
        allPredOutputs(iCL,:) = shiftOutputs*ones(1,length(nComps));
        if useVIP
            allPredOutputs_vip(iCL,:,:,:) = shiftOutputs*ones(length(nComps),nVIP,length(nComps2),'single');
            nVar_vip(iCL,:,:,:) = 0;
        end
        continue
    end

    [~,YL,XS,~,~,W] = plsregress_simpls(NormInput, NormOutputs, max(nComps));

    R2 = zeros(max(nComps),1);
    nWei = zeros(1,max(nComps));
    if doVIP
        for j=1:max(nComps); nWei(j) = norm(W(:,j)); end
        for i=1:nVar
            if i==1;
                for j=1:max(nComps)
                    R2(j) = sum( (YL(:,j).^2) * ...
                        (XS(:,j)'*XS(:,j)) ) ;
                end;end
        end
    end

    tnVar_vip = zeros(length(nComps),nVIP,length(nComps2),'single');
    tallPredOutputs_vip = zeros(length(nComps),nVIP,length(nComps2),'single');
    tallbeta = zeros(size(Inputs,2),length(nComps));
    tallVIP = tallbeta;
    tallPredOutputs = NaN(1,length(nComps));


    for j=1:length(nComps)  % prediction for the different number of components
        n = nComps(j);

        beta = W(:,1:n)*YL(:,1:n)';
        tallbeta(:,j) = beta;
        tallPredOutputs(j) = ( Inputs(iCL,:) -shiftInputs)*beta + shiftOutputs;

        tVIP = zeros(nVar,1);
        %%%%%%% runing the optimization of the Inputs using the VIP
        if doVIP

            % calculating the VIP for the training set
            for i=1:nVar
                tVIP(i) = sqrt( nVar* (((W(i,1:n)./nWei(1:n)).^2)*R2(1:n)) /sum(R2(1:n)) );
            end
            % recording the VIP
            tallVIP(:,j) = tVIP;

        end % end of the VIP

        if useVIP

            for i = 1:nVIP
                VIPselected = (tVIP >= allVIPthres(i));
                if ~any(VIPselected)
                    break
                end
                tnComp = max(1,min(max(nComps2),(rank(NormInput(:,VIPselected))-3)));
                [~,YL2,~,~,~,W2] = plsregress_simpls(NormInput(:,VIPselected), NormOutputs, tnComp);

                for n2=intersect(1:tnComp,nComps2)  % prediction for the different number of components
                    k = find(nComps2==n2);
                    beta = W2(:,1:n2)*YL2(:,1:n2)';
                    tnVar_vip(n,i,k) = sum(VIPselected);
                    tallPredOutputs_vip(n,i,k) = ...
                        ( Inputs(iCL,VIPselected) -shiftInputs(VIPselected)) ...
                        *beta + shiftOutputs;
                end
            end

        end
    end

    allPredOutputs(iCL,:) = tallPredOutputs;
    if useVIP
        allPredOutputs_vip(iCL,:,:,:) = tallPredOutputs_vip;
        nVar_vip(iCL,:,:,:) = tnVar_vip;
    end

end

allMSE = zeros(1,length(nComps));
for n=1:length(nComps)
    allMSE(n) = sum((allPredOutputs(:,n) - Outputs).^2);
end
details.allMSE = allMSE;

details.optnComp = find(allMSE <= (min(allMSE)+std(allMSE)),1,'first');

SSquares = sum( (Outputs-mean(Outputs)).^2 );

chosen_q2 = 1 - allMSE(details.optnComp)/SSquares;

NormInput = Inputs-repmat(mean(Inputs),nCL,1);
NormOutputs = Outputs-mean(Outputs);
[~,YL,XS,~,~,W] = plsregress_simpls(NormInput, NormOutputs, details.optnComp);

B = W*YL';
fit_output = NormInput*B + mean(Outputs);
r2 = 1 -  sum( (fit_output-Outputs).^2 ) /SSquares;

if doVIP
    optnComp = details.optnComp;

    % generate the full model without leave-1-out
    R2 = zeros(optnComp,1);
    nWei = zeros(1,optnComp);
    if doVIP
        for j=1:optnComp; nWei(j) = norm(W(:,j)); end
        for i=1:nVar
            if i==1;
                for j=1:optnComp
                    R2(j) = sum( (YL(:,j).^2) * ...
                        (XS(:,j)'*XS(:,j)) ) ;
                end;end
        end
    end
    details.VIP = zeros(nVar,1);
    % calculating the VIP for the training set
    for i=1:nVar
        details.VIP(i) = sqrt( nVar* (((W(i,:)./nWei).^2)*R2) /sum(R2) );
    end
end

if useVIP
    allMSE_vip = NaN(length(nComps),nVIP,length(nComps2));
    for n=1:length(nComps)
        for i = 1:nVIP
            for n2=1:length(nComps2)
                if all(isnan(allPredOutputs_vip(:,n,i,n2)))
                    break
                end
                allMSE_vip(n,i,n2) = sum((allPredOutputs_vip(:,n,i,n2) - Outputs).^2);
            end
        end
    end


    % chossign the parameter as minimal # of VIP then minimal number of
    % Components
    VIPidx = [];
    for i = nVIP:-1:1
        [nidx n2idx] = find(allMSE_vip(:,i,:)<=quantile(reshape(allMSE_vip,1,[]),.01));
        for j=1:length(nidx)
            if ~isempty(nidx) && min(nVar_vip(:,nidx(j),i,n2idx(j)))>0
                nidx = nidx(j);
                VIPidx = i;
                n2idx = n2idx(j);
                break
            end
        end
        if ~isempty(VIPidx)
            break
        end
    end
    if isempty(VIPidx)
        for i = nVIP:-1:1
            [nidx n2idx] = find(allMSE_vip(:,i,:)<=(1e-10+quantile( ...
                reshape(allMSE_vip+ shiftdim( (1+max(allMSE_vip(:)))*(min(nVar_vip,[],1)==0) ,1) ,1,[]) ,.01) ));
            for j=1:length(nidx)
                if ~isempty(nidx) && min(nVar_vip(:,nidx(j),i,n2idx(j)))>0
                    nidx = nidx(j);
                    VIPidx = i;
                    n2idx = n2idx(j);
                    break
                end
            end
            if ~isempty(VIPidx)
                break
            end
        end
    end
    if isempty(VIPidx)
        VIPidx = 1;
        nidx = min(size(allMSE_vip,1),2);
        n2idx = min(size(allMSE_vip,3),2);
    end


    details_vip.nVar_vip = nVar_vip;
    details_vip.allMSE_vip = allMSE_vip;
    details_vip.optnComp = nComps(nidx);
    details_vip.VIPcutoff = allVIPthres(VIPidx);
    details_vip.optnComp2 = nComps2(n2idx);

    chosen_q2_vip = 1 - allMSE_vip(nidx, VIPidx, n2idx) / SSquares;

    % generating the final model
    % generate the full model without leave-1-out
    [~,YL,XS,~,~,W] = plsregress_simpls(NormInput, NormOutputs, details_vip.optnComp);
    R2 = zeros(details_vip.optnComp,1);
    nWei = zeros(1,details_vip.optnComp);
    if doVIP
        for j=1:details_vip.optnComp; nWei(j) = norm(W(:,j)); end
        for i=1:nVar
            if i==1;
                for j=1:details_vip.optnComp
                    R2(j) = sum( (YL(:,j).^2) * ...
                        (XS(:,j)'*XS(:,j)) ) ;
                end;end
        end
    end
    details_vip.VIP = zeros(nVar,1);
    % calculating the VIP for the training set
    for i=1:nVar
        details_vip.VIP(i) = sqrt( nVar* (((W(i,:)./nWei).^2)*R2) /sum(R2) );
    end

    Varselected = (details_vip.VIP >= details_vip.VIPcutoff);
    if sum(Varselected)>0

        if details_vip.optnComp2>=rank(NormInput(:,Varselected))
            details_vip.optnComp2 = max(1,rank(NormInput(:,Varselected))-1);
        end
        try
            [~,YL,~,~,~,W] = plsregress_simpls(NormInput(:,Varselected), NormOutputs, details_vip.optnComp2);
        catch
            if details_vip.optnComp2>1
                disp(['Error; changeing to ' num2str(details_vip.optnComp2-1) ' components'])
                [~,YL,~,~,~,W] = plsregress_simpls(NormInput(:,Varselected), NormOutputs, details_vip.optnComp2-1);
            else
                disp(['Error; changeing to ' num2str(details_vip.optnComp2-1) ' components -> set to ones.'])
                YL = ones(1,1);
                W = ones(sum(Varselected),1);
            end
        end
        B_vip = 0*B;
        B_vip(Varselected) = W*YL';
        fit_output_vip = NormInput(:,Varselected)*(W*YL') + mean(Outputs);

    else
        B_vip = 0*B;
        fit_output_vip = mean(Outputs);
    end

    r2_vip = 1 -  sum( (fit_output_vip-Outputs).^2 ) /SSquares;
end
