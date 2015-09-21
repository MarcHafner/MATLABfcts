plate_design = Read_Plate_Design('example_plateDesign.tsv');

%%
Raw_data = Load_Plate_Data(plate_design);

%%
MOMP_data = Evaluate_MOMPTime(Raw_data);

%%
Process_data = Process_trajectories(MOMP_data);

%%
Fitted_data = FitQuadratic(Process_data);
