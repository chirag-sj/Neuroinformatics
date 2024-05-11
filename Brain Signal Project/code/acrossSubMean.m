function [fz_1,fz_2,fz_3,fcz_1,fcz_2,fcz_3,cz_1,cz_2,cz_3] = acrossSubMean(subj_id_order,final_data,grp)
%ACROSSSUBMEAN Summary of this function goes here
%   Detailed explanation goes here
fz_1=[];
for i=1:length(subj_id_order)
    if ismember(subj_id_order(i),grp)
        strct_temp = final_data(subj_id_order(i));
        if isempty(fz_1)
            fz_1 = strct_temp{1,1}(1).FC3;
            fz_2 = strct_temp{1,1}(2).FC3;
            fz_3 = strct_temp{1,1}(3).FC3;
            
            fcz_1 = strct_temp{1,1}(1).FC4;
            fcz_2 = strct_temp{1,1}(2).FC4;
            fcz_3 = strct_temp{1,1}(3).FC4;
            
            cz_1 = strct_temp{1,1}(1).C3;
            cz_2 = strct_temp{1,1}(2).C3;
            cz_3 = strct_temp{1,1}(3).C3;

            disp(isequal(fz_1,fcz_1))
        else
            fz_1 = (fz_1+strct_temp{1,1}(1).FC3);
            fz_2 = (fz_2+strct_temp{1,1}(2).FC3);
            fz_3 = (fz_3+strct_temp{1,1}(3).FC3);
        
            fcz_1 = (fcz_1+strct_temp{1,1}(1).FC4);
            fcz_2 = (fcz_2+strct_temp{1,1}(2).FC4);
            fcz_3 = (fcz_3+strct_temp{1,1}(3).FC4);
        
            cz_1 = (cz_1+strct_temp{1,1}(1).C3)./2;
            cz_2 = (cz_2+strct_temp{1,1}(2).C3)./2;
            cz_3 = (cz_3+strct_temp{1,1}(3).C3)./2;
        end
    end
end
end

