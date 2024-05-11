%The ERPs are calculated by averaging across trials for every sample in the time series, separately for each subject, electrode, and condition.
% subject_ind = 1;
% electrode = 'FCz';
% condition_val = 2;
% % [sub_data,demographic_sub,ERPdata_sub,mergedTrialData_sub,time_sub] = loadData(subject_ind);
% 
% trial_data = [];
% for i=1:100
%     if isempty(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).(electrode))
%         continue
%     elseif isempty(trial_data)
%         trial_data = sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).(electrode)';
%     else
%         sub_trial_data = sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).(electrode)';
%         trial_data = cat(1, trial_data,sub_trial_data);
%     end
% end
% 
% sub_elec_cond_erp = mean(trial_data,1);


% HCsubject_ind = 1;
% columnLabels = readtable("../archive/columnLabels.csv").Properties.VariableNames;
% HCpath = strcat("../archive/",num2str(HCsubject_ind),".csv/",num2str(HCsubject_ind),".csv");
% HCsub_data = readtable(HCpath);
% HCsub_data.Properties.VariableNames = columnLabels;
% 
% SZsubject_ind = 67;
% SZpath = strcat("../archive/",num2str(SZsubject_ind),".csv/",num2str(SZsubject_ind),".csv");
% SZsub_data = readtable(SZpath);
% SZcolumnLabels = readtable("../archive/columnLabels.csv").Properties.VariableNames;
% SZsub_data.Properties.VariableNames = columnLabels;


[fz_erp_HC_c1,fcz_erp_HC_c1,cz_erp_HC_c1] = getERP(HCsub_data,1);
[fz_erp_HC_c2,fcz_erp_HC_c2,cz_erp_HC_c2] = getERP(HCsub_data,2);
[fz_erp_HC_c3,fcz_erp_HC_c3,cz_erp_HC_c3] = getERP(HCsub_data,3);

[fz_erp_SZ_c1,fcz_erp_SZ_c1,cz_erp_SZ_c1] = getERP(SZsub_data,1);
[fz_erp_SZ_c2,fcz_erp_SZ_c2,cz_erp_SZ_c2] = getERP(SZsub_data,2);
[fz_erp_SZ_c3,fcz_erp_SZ_c3,cz_erp_SZ_c3] = getERP(SZsub_data,3);

time_axis = readmatrix('../archive/time.csv');
time_axis = time_axis(:,2);

figure
subplot(121)
plot(time_axis,fz_erp_HC_c1,'b')
hold on
plot(time_axis,fz_erp_HC_c2,'m') 
% hold on
% plot(time_axis,fz_erp_HC_c3,'b')
title('Healthy Control')

subplot(122)
plot(time_axis,fz_erp_SZ_c1,'b')
hold on
plot(time_axis,fz_erp_SZ_c2,'m')
% hold on
% plot(time_axis,fz_erp_SZ_c3,'b')
title('Schizophrenic')

figure
subplot(121)
plot(time_axis,fcz_erp_HC_c1,'b')
hold on
plot(time_axis,fcz_erp_HC_c2,'m') 
% hold on
% plot(time_axis,fcz_erp_HC_c3,'b')
title('Healthy Control')

subplot(122)
plot(time_axis,fcz_erp_SZ_c1,'b')
hold on
plot(time_axis,fcz_erp_SZ_c2,'m')
% hold on
% plot(time_axis,fcz_erp_SZ_c3,'b')
title('Schizophrenic')
legend('Button Tone','Play Tone');
figure
subplot(121)
plot(time_axis,cz_erp_HC_c1,'b')
hold on
plot(time_axis,cz_erp_HC_c2,'m') 
% hold on
% plot(time_axis,cz_erp_HC_c3,'b')
title('Healthy Control')

subplot(122)
plot(time_axis,cz_erp_SZ_c1,'b')
hold on
plot(time_axis,cz_erp_SZ_c2,'m')
% hold on
% plot(time_axis,cz_erp_SZ_c3,'b')
title('Schizophrenic')

legend('Button Tone','Play Tone');