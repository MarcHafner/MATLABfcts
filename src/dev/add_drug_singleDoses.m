function drugs_struct_out = add_drug_singleDoses(drugs_struct, new_drug_struct)

assert(~ismember(new_drug_struct.name, {drugs_struct.name}))

new_drug_struct.volume = new_drug_struct.well_volume*sum(new_drug_struct.SingleDoses)/...
        (1e3*new_drug_struct.nominal_conc);
new_drug_struct.Doses = [];
new_drug_struct.layout = zeros(16,24);


drugs = NaN(16,24,length(drugs_struct));
for i=1:length(drugs_struct)
    drugs(:,:,i) = drugs_struct(i).layout;
end

current_ctrl = all(drugs==0,3);
%%
ctrl_idx = find(current_ctrl);
ctrl_idx = ctrl_idx(randperm(length(ctrl_idx)));

%%

for i=1:length(new_drug_struct.SingleDoses)
    new_drug_struct.layout(ctrl_idx(i)) = new_drug_struct.SingleDoses(i);
end


drugs_struct_out = drugs_struct;
drugs_struct_out(end+1) = new_drug_struct;
