function drugs_struct_out = force_ctrl_pos(drugs_struct, ctrl_wells)

drugs_struct_out = drugs_struct;

drugs = NaN(16,24,length(drugs_struct));
for i=1:length(drugs_struct)
    drugs(:,:,i) = drugs_struct(i).layout;
end

current_ctrl = all(drugs==0,3);
%%
ctrl_idx = find(current_ctrl);

forced_ctrls = NaN(length(ctrl_wells),2);
ctrl_wells = upper(ctrl_wells);
for i=1:length(ctrl_wells)
    forced_ctrls(i,1) = ctrl_wells{i}(1)-64;
    forced_ctrls(i,2) = str2num(ctrl_wells{i}(2:end));
end

forced_ctrls = sub2ind([16 24], forced_ctrls(:,1), forced_ctrls(:,2));

inter_idx = intersect(forced_ctrls,ctrl_idx);
forced_ctrls = setdiff(forced_ctrls,inter_idx);
ctrl_idx = setdiff(ctrl_idx,inter_idx);

% random swap
ctrl_idx = ctrl_idx (randperm(length(ctrl_idx))<=length(forced_ctrls));



for i=1:length(drugs_struct)
    for j=1:length(ctrl_idx)
        drugs_struct_out(i).layout(ctrl_idx(j)) = drugs_struct(i).layout(forced_ctrls(j));
        drugs_struct_out(i).layout(forced_ctrls(j)) = 0;
    end
end
