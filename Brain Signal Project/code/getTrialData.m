function [trial_data] = getTrialData(data,trial_ind)
%GETTRIALDATA Summary of this function goes here
%   Detailed explanation goes here
trial_data = data(data.trial==trial_ind,:);
end

