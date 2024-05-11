listCSV = dir('../*/*/*.csv');
listCSV_folder = {listCSV.folder};
listCSV_name = {listCSV.name};
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

% listInfoCSV = dir('../archive/*.csv');
% listFile = listInfoCSV([listInfoCSV.bytes]>0);
% listFile_folder = {listFile.folder};
% listFile_name = {listFile.name};
% load(append(listFile_folder(1),listFile_name(1)))

for i=1:length(subj_id_order)
    path2csv = append('../archive/',int2str(subj_id_order(i)),'.csv/',int2str(subj_id_order(i)),'.csv');
end
final_data_new = {};

for i=1:length(subj_id_order)
    disp(i)
    path = strcat("../archive/",num2str(subj_id_order(i)),".csv/",num2str(subj_id_order(i)),".csv");
    sub = readtable(path);
    columnLabels = readtable("../archive/columnLabels.csv").Properties.VariableNames;
    sub.Properties.VariableNames = columnLabels;
    subj_data = {};
    for j=1:3
        [fz,fcz,cz] = getERP(sub,j);
        subj_data(j).('FC3') = fz;
        subj_data(j).('FC4') = fcz;
        subj_data(j).('C3') = cz;
    end
    final_data_new{subj_id_order(i)} = subj_data;
    
end

