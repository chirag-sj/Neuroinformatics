load final_data.mat
load all_TF.mat


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


[~,indHC,~] = intersect(subj_id_order,grpHC);
[~,indSZ,~] = intersect(subj_id_order,grpSz);
% Define the number of permutations
num_permutations = 1000;

% Initialize arrays to store observed and permuted test statistics
observed_stat = zeros(size(storeAllTF, 2), size(storeAllTF, 3)); % Assuming storeAllTF is subject x frequency x time
permuted_stat = zeros(num_permutations, size(storeAllTF, 2), size(storeAllTF, 3));

% Compute observed test statistic
% For example, you can calculate the difference in means between the two groups
observed_diff = mean(storeAllTF(indSZ, :, :), 1) - mean(storeAllTF(indHC, :, :), 1);

% Perform permutations
for perm = 1:num_permutations
    % Shuffle group labels
    shuffled_labels = datasample([ones(length(indSZ), 1); 2*ones(length(indHC), 1)], length(indSZ) + length(indHC), 'Replace', false);
    
    % Compute test statistic for permuted data
    permuted_diff = mean(storeAllTF(shuffled_labels == 1, :, :), 1) - mean(storeAllTF(shuffled_labels == 2, :, :), 1);
    permuted_stat(perm, :, :) = permuted_diff;
end

% Calculate p-values for each frequency-time point
p_values = zeros(size(observed_diff)); % Preallocate p_values matrix
for freq = 1:size(observed_diff, 1)
    for time = 1:size(observed_diff, 2)
        % Perform t-test for each frequency-time point pair
        [~, p_values(freq, time)] = ttest(squeeze(observed_diff(freq, time, :)));
    end
end

% Plot p-values
figure;
imagesc(squeeze(p_values));
colorbar; % Add color bar
title('p-values'); % Add title
xlabel('Time'); % Add x-axis label
ylabel('Frequency'); % Add y-axis label

