%% Load data.
sub1 = readtable("../archive/1.csv/1.csv");
columnLabels = readtable("../archive/columnLabels.csv").Properties.VariableNames;
sub1.Properties.VariableNames = columnLabels;
demographic= readtable("../archive/demographic.csv");
ERPdata = readtable("../archive/ERPdata.csv");
mergedTrialData = readtable("../archive/mergedTrialData.csv");
time = readtable("../archive/time.csv");
% Important GetInfo
%unique(sub1.condition(find(sub1.trial==trial_ind)))
%length(find(sub1.condition==cond_ind))