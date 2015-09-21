% load and analyze the single cell data
%
%

function output = LoadFitData(filename, savefilename)

plate_design = Read_Plate_Design(filename);

if ~exist('savefilename','var')
    savefilename = ['./datasets/' plate_design.tag '.mat'];
end

Raw_data = Load_Plate_Data(plate_design);
%

MOMP_data = Evaluate_MOMPTime(Raw_data);
save(savefilename, 'MOMP_data')

Process_data = Process_trajectories(MOMP_data);
Fitted_data = FitQuadratic_Ftest(Process_data);
save(savefilename, 'Fitted_data')

if exist(savefilename, 'file')
    output = 0;
else
    output = -1;
end
