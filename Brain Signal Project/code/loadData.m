function [sub,demographic,ERPdata,mergedTrialData,time] = loadData(subject_ind)
%LOADDATA Summary of this function goes here
%   Detailed explanation goes here
path = strcat("../archive/",num2str(subject_ind),".csv/",num2str(subject_ind),".csv");
sub = readtable(path);
columnLabels = readtable("../archive/columnLabels.csv").Properties.VariableNames;
sub.Properties.VariableNames = columnLabels;
demographic= readtable("../archive/demographic.csv");
demographic = demographic(demographic.subject==1,:);
ERPdata = readtable("../archive/ERPdata.csv");
ERPdata = ERPdata(ERPdata.subject==subject_ind,:);
mergedTrialData = readtable("../archive/mergedTrialData.csv");
mergedTrialData = mergedTrialData(mergedTrialData.subject==subject_ind,:);
time = readtable("../archive/time.csv");
end