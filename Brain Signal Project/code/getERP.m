function [erpfz,erpfcz,erpcz] = getERP(sub_data ,condition_val )

%GETERP Summary of this function goes here
% Detailed explanation goes here
% The ERPs are calculated by averaging across trials for every sample in the time series, separately for each subject, electrode, and condition.
% subject_ind = 1;
% electrode = 'FCz';
% condition_val = 2;
% [sub_data,demographic_sub,ERPdata_sub,mergedTrialData_sub,time_sub] = loadData(subject_ind);
% srate = 1024;
% nyquist          = srate/2;
% transition_width = 0.15; % percent
% % finally, filter from 5-15
% filter_low    =  0.5; % Hz
% filter_high   = 15; % Hz
% ffrequencies  = [ 0 filter_low*(1-transition_width) filter_low filter_high filter_high*(1+transition_width) nyquist ]/nyquist;
% idealresponse = [ 0 0 1 1 0 0 ];
% filterweights = firls(round(3*(srate/filter_low)),ffrequencies,idealresponse);
% erp_5to15     = filtfilt(filterweights,1,double(erp));



trial_data_Fc = [];
trial_data_FCz = [];
trial_data_Cz = [];


for i=1:100
    if isempty(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).('FC3'))
        continue
    elseif isempty(trial_data_Fc)
        trial_data_Fc = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).('FC3')');
        trial_data_FCz = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).('FC4')');
        trial_data_Cz = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).('C3')');
    else
        sub_trial_data_Fc = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).('FC3')');
        trial_data_Fc = cat(1, trial_data_Fc,sub_trial_data_Fc);

        sub_trial_data_FCz = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).('FC4')');
        trial_data_FCz = cat(1, trial_data_FCz,sub_trial_data_FCz);

        sub_trial_data_Cz = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).('C3')');
        trial_data_Cz = cat(1, trial_data_Cz,sub_trial_data_Cz);
    end
end

% sub_elec_cond_erp = mean(trial_data_Fc,1);
% erpfz = filtfilt(filterweights,1,double(mean(trial_data_Fc,1)));
% erpfcz = filtfilt(filterweights,1,double(mean(trial_data_FCz,1)));
% erpcz = filtfilt(filterweights,1,double(mean(trial_data_Cz,1)));
erpfz = mean(trial_data_Fc,1);
erpfcz = mean(trial_data_FCz,1);
erpcz = mean(trial_data_Cz,1);

% fs = 1024;
% f_low = 0.5;
% f_high = 15;
% [b, a] = butter(2, [f_low/(fs/2), f_high/(fs/2)], 'bandpass'); 
% 
% erpfz = filtfilt(b, a, erpfz);
% erpfcz = filtfilt(b, a, erpfcz);
% erpcz = filtfilt(b, a, erpcz);

% figure
% plot(sub_elec_cond_erp)
end

