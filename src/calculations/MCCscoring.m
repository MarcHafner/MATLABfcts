function score = MCCscoring (ConfMx)
% score = MCCscoring (ConfMx)
%   ConfMx =
%       [  TP    FP
%          FN    TN ];
%

TP = ConfMx(1,1);
TN = ConfMx(2,2);
FN = ConfMx(2,1);
FP = ConfMx(1,2);


if (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)~=0
    score = ( (TP*TN)-(FP*FN) )/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
else
    score = (TP*TN)-(FP*FN);
end

if score<-1 || score>1
    disp([score TP TN FN FP])
    error('mcc out of range')
end
