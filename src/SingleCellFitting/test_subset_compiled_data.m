% eample final processing

final_data = filter_DeltaMOMP(compiled_data);

% example for choosing some conditions
Idx = find(final_data.ExpKey.TRAIL==25 & final_data.ExpKey.ABT==0 );

% extracting the specific conditions
Exp_data = subset_compiled_data(final_data, Idx);

% removing the trajectories after MOMP
Exp_data = remove_postMOMP(Exp_data);
