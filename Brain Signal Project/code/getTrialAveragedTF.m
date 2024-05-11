
load final_data.mat

% [~,indHC,~] = intersect(subj_id_order,grpHC);
min_freq = 2;
max_freq = 64;
num_frex = 30;

% other wavelet parameters
srate=1024;
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);

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
storeAllTF = zeros(40,30,length(timeaxis));
for i=1:length(subj_id_order)
    strct_temp = final_data(subj_id_order(i));
    avgSubCndElec = strct_temp{1,1}(1).fcz;
    subjTF = getTF(avgSubCndElec);
    storeAllTF(i,:,:)=subjTF;
end

% Define plot parameters
num_rows = 10; % Number of rows in the subplot grid
num_cols = 4;  % Number of columns in the subplot grid

% % Get the number of time points
% num_time_points = length(timeaxis);
% ytickskip = round(logspace(log10(1), log10(length(frequencies)), 10));
% yticklabels = round(frequencies(ytickskip));
% % Loop through the images and plot them
% figure;
% for i = 1:size(storeAllTF, 1)
%     subplot(num_rows, num_cols, i);
%     imagesc(timeaxis, frequencies, squeeze(storeAllTF(i, :, :))); % Display the image with timeaxis as x-axis
%     colormap jet; % Set the colormap
%     colorbar; % Show colorbar
%     
%     % Determine subject index and group
%     subj_idx = subj_id_order(i);
%     if ismember(subj_idx, grpSz)
%         group_label = 'Schizophrenic';
%     elseif ismember(subj_idx, grpHC)
%         group_label = 'Healthy Control';
%     else
%         group_label = 'Unknown Group';
%     end
%     title(['Subject ', num2str(subj_idx), ': ', group_label]); % Set subplot title
%     xlabel('Time (ms)'); % Set x-axis label
%     ylabel('Frequency'); % Set y-axis label
%     xlim([-100 300]);
%     set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-100 300],'clim',[-12 12])
% end
% 
% % Adjust the layout
% sgtitle('Images from storeAllTF Array'); % Set the main title


num_rows = 10; 
num_cols = 4; 

num_time_points = length(timeaxis);

ytickskip = round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100;

figure;
for i = 1:size(storeAllTF, 1)
    subplot(num_rows, num_cols, i);
    imagesc(timeaxis, frequencies, squeeze(storeAllTF(i, :, :))); 
    colormap jet; 
    colorbar; 

    subj_idx = subj_id_order(i);
    if ismember(subj_idx, grpSz)
        group_label = 'Schizophrenic';
    elseif ismember(subj_idx, grpHC)
        group_label = 'Healthy Control';
    else
        group_label = 'Unknown Group';
    end
    title(['Subject ', num2str(subj_idx), ': ', group_label]); 
    xlabel('Time (ms)'); 
    ylabel('Frequency');
    xlim([-500 500]); 
    set(gca, 'ytick', round(ytickskip),  'ydir', 'normal', 'clim', [-12 12]); 
end
% 'yticklabel', yticklabels,
sgtitle('Time Frequency Analysis'); 
