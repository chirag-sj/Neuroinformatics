A = rand(4,8);
mat = A_func(A);

columnLabels = {'Row Index', 'Column Index', 'Label (0<0.5<1)'};  % Replace with your actual labels
dataTable = array2table(mat, 'VariableNames', columnLabels);
writetable(dataTable, 'output_file.txt', 'Delimiter', '\t');
