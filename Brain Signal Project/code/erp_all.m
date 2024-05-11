load("final_data_new.mat");

listCSV = dir('../*/*/*.csv');
listCSV_folder = {listCSV.folder};
listCSV_name = {listCSV.name};
timeaxis = readtable("../archive/time.csv");
timeaxis = timeaxis.time_ms;
subj_path_all = {};
subj_id_order = [];
for i = 1:length(listCSV)
    temp_listCSV_folder = listCSV_folder(i);
    temp_listCSV_name = listCSV_name(i);
    subj_path = append(temp_listCSV_folder,temp_listCSV_name);
    subj_path_all(i)=subj_path;
    split_temp_indx = split(temp_listCSV_name,'.');
    split_temp_indx = split_temp_indx(1);
    subj_id_order(i) = str2num(split_temp_indx{1});
end

demographic= readtable("../archive/demographic.csv");

grpHC = find(demographic.group==0);
grpSz = find(demographic.group==1);

[fz_1_hc,fz_2_hc,fz_3_hc,fcz_1_hc,fcz_2_hc,fcz_3_hc,cz_1_hc,cz_2_hc,cz_3_hc] = acrossSubMean(subj_id_order, final_data_new,grpHC);
[fz_1_sz,fz_2_sz,fz_3_sz,fcz_1_sz,fcz_2_sz,fcz_3_sz,cz_1_sz,cz_2_sz,cz_3_sz] = acrossSubMean(subj_id_order,final_data_new,grpSz);



figure

subplot(1,2,1)
plot(timeaxis,fz_1_hc,'b',timeaxis,fz_2_hc,'m',timeaxis,fz_3_hc,'k')
title('Healthy Control')
set(gca,'xlim',[-100 300])
subplot(1,2,2)
plot(timeaxis,fz_1_sz,'b',timeaxis,fz_2_sz,'m',timeaxis,fz_3_sz,'k')
title('Schizophrenic')
set(gca,'xlim',[-100 300])
legend('Button Tone','Play Tone','Button Alone')
figure
subplot(1,2,1)
plot(timeaxis,fcz_1_hc,'b',timeaxis,fcz_2_hc,'m',timeaxis,fcz_3_hc,'k')
title('Healthy Control')
set(gca,'xlim',[-100 300])
subplot(1,2,2)
plot(timeaxis,fcz_1_sz,'b',timeaxis,fcz_2_sz,'m',timeaxis,fcz_3_sz,'k')
title('Schizophrenic')
set(gca,'xlim',[-100 300])
legend('Button Tone','Play Tone','Button Alone')
figure
subplot(1,2,1)
plot(timeaxis,cz_1_hc,'b',timeaxis,cz_2_hc,'m',timeaxis,cz_3_hc,'k')
title('Healthy Control')
set(gca,'xlim',[-100 300])
subplot(1,2,2)
plot(timeaxis,cz_1_sz,'b',timeaxis,cz_2_sz,'m',timeaxis,cz_3_sz,'k')
title('Schizophrenic')
set(gca,'xlim',[-100 300])
legend('Button Tone','Play Tone','Button Alone')