function drugs_struct = ImportMap_xlsx(filename, savefilename, sheetname)


if exist('sheetname','var')
    [~,txt,data] = xlsread(filename,sheetname);
else
    [~,txt,data] = xlsread(filename);
end

%%
Drugs = unique(setdiff(txt(:,1),[num2cell('A':'P') {''}]));


drugs_struct = struct('name', Drugs, 'nominal_conc', NaN(size(Drugs)), ...
    'layout',zeros(16,24));

alllayouts = zeros(16,24,length(Drugs));

for iD=1:length(Drugs)
    Didx = find(strcmp(txt(:,1), Drugs{iD}));

    assert(data{Didx+2,1}=='A')
    assert(data{Didx+17,1}=='P')
    assert(data{Didx+1,2}==1)
    assert(data{Didx+1,25}==24)

    drugs_struct(iD).nominal_conc = data{Didx,2};
    drugs_struct(iD).layout = cell2mat(data(Didx+(2:17), 2:25));

    alllayouts(:,:,iD) = cell2mat(data(Didx+(2:17), 2:25));
end

%%

tested_pairs = zeros(0,2);

for iD=1:length(Drugs)

    NoOtherDrug = all(alllayouts(:,:,setdiff(1:length(Drugs),iD))==0,3);
    drugs_struct(iD).Doses = unique(drugs_struct(iD).layout(~NoOtherDrug(:)))';
    drugs_struct(iD).SingleDoses = unique(drugs_struct(iD).layout(NoOtherDrug(:)))';

    for iD2=(iD+1):length(Drugs)
        if any(any(alllayouts(:,:,iD)>0 & alllayouts(:,:,iD2)>0))
            tested_pairs = [tested_pairs; iD iD2];
        end
    end
end

n = hist(tested_pairs(:),1:length(Drugs));

Complete = false;
if length(unique(n))>2
    disp('No pattern with primary vs secondary drugs')
elseif length(unique(n))==2
    Primary = n==max(n);
    Complete = true;
    for i=find(Primary)
        temp = tested_pairs( any(tested_pairs==i,2),:);
        temp = setdiff(temp(:),i)';
        if length(temp)~=sum(~Primary) || any(temp~=find(~Primary))
            fprintf('%s is not tested against all and only secondary Drugs\n',...
                Drugs{i});
            Complete = false;
            break
        end
    end
end

if Complete
    for iD=1:length(Drugs)
        drugs_struct(iD).Primary = Primary(iD);
    end
    drugs_struct = drugs_struct(sortidx([drugs_struct.Primary],'descend'));
end

if exist('savefilename','var')
    if isempty(savefilename)
        savefilename = filename(1:find(filename=='.',1,'last'));
        savefilename = [savefilename 'mat'];
    end
    save(savefilename, 'drugs_struct')
end
