% This code evaluates the Deep ensemble transfer learning (DETL) model output and
% fuses results from CC and MLO views of the validation data.
% Deep learning model output is taken from the jupyter notebook file.
% val1 folder will be used for validation

% res_table: first column includes the combined estimation of DETL from CC
% and/or MLO views
% second column is 0 if the leison is benign and 1 if the lesion is malign

clear; close all; clc; format compact

% Read BCDR table in excel format
bcdr_table = readtable('bcdr_precomputed_features_filenames_30mar22.xlsx');

% Load train validation test data indices. This is for consistency
% between experiments.
load bcdr_CC_fuse_MLO_5f_xval_tst

%Load validation patch names and the corresponding deep ensemble model
%predictions
load testPatches_15feb23.mat 


% Here the patient ids are given in the order of the bcdr_rmid_table
alphacc=0.55; %The weight coefficient for CC view. I expect alphacc>0.5 
% would give better results than equal contribution
rc=0;

% Use tst_set 
for i1=1:length(tst_set)
    %Find actual lesion id
    les_id=bcdr_rmid_table(tst_set(i1),3);
    %Find row numbers of this lesion in the BCDR table
    les_rows=find(table2array(bcdr_table(:,3))==les_id);
    if length(les_rows)==1
        rc=rc+1;
        %convert the patch name
        bcdr_overlay_name=char(string(table2array(bcdr_table(les_rows(1),20))));
        bcdr_overlay_name=strcat(bcdr_overlay_name(1,1:end-4),'_patch.png');
        %find patch name row number in the prediction column
        row_num=0;
        for j1=1:length(patchName)
            pred_patch_name=char(patchName(j1));
            if strcmp(pred_patch_name(6:end),bcdr_overlay_name)
                row_num=j1;
                pred_patch_name=char(patchName(j1));
                break
            end
        end
        if row_num==0
            error('patch not found');
        end
        %read deep ensemble model prediction and write to the result table
        %Hint: Here we do not need to combine results 
        res_table(rc,1)=Y_pred_XNR(row_num,1);
        % write the ground truth
        if strcmp(pred_patch_name(1:3),'mal')
            res_table(rc,2)=1;
        end
        if strcmp(pred_patch_name(1:3),'ben')
            res_table(rc,2)=0;
        end
        clear pred_patch_name
        
    %elseif length(les_rows)==2   
    else %combine CC and MLO results. Here there is always two views
        rc=rc+1;
        %Read first overlay from BCDR table
        overlay_name1=char(string(table2array(bcdr_table(les_rows(1),20))));
        %find view type
        view_type1=char(string(table2array(bcdr_table(les_rows(1),6))));
        view_type1=view_type1(end-5:end);
        %Read second overlay from BCDR table
        overlay_name2=char(string(table2array(bcdr_table(les_rows(2),20))));
        view_type2=char(string(table2array(bcdr_table(les_rows(2),6))));
        view_type2=view_type2(end-5:end);
        
        % Read prediction1 for patch1 from DETL
        bcdr_overlay_name1=strcat(overlay_name1(1,1:end-4),'_patch.png');
        %find patch name row number in the prediction column
        row_num=0;
        for j1=1:length(patchName)
            pred_patch_name=char(patchName(j1));
            if strcmp(pred_patch_name(6:end),bcdr_overlay_name1)
                row_num=j1;
                pred_patch_name=char(patchName(j1));
                break
            end
        end
        %Check for error
        if row_num==0
           error('patch not found');
        end
        res1=Y_pred_XNR(row_num,1);
        
        % Read prediction2 for patch2 from DETL
        bcdr_overlay_name2=strcat(overlay_name2(1,1:end-4),'_patch.png');
        %find patch name row number in the prediction column
        row_num=0;
        for j1=1:length(patchName)
            pred_patch_name=char(patchName(j1));
            if strcmp(pred_patch_name(6:end),bcdr_overlay_name2)
                row_num=j1;
                pred_patch_name=char(patchName(j1));
                break
            end
        end
        %Check for error
        if row_num==0
           error('patch not found');
        end
        res2=Y_pred_XNR(row_num,1);
       
        % Calculate result by combining CC and MLO
        if strcmp(view_type1,'caudal')
            res=res1*alphacc+res2*(1-alphacc);
        else
            res=res1*(1-alphacc)+res2*alphacc;
        end
        
        % Write to the result table
        res_table(rc,1)=res;
        % write the ground truth
        if strcmp(pred_patch_name(1:3),'mal')
            res_table(rc,2)=1;
        end
        if strcmp(pred_patch_name(1:3),'ben')
            res_table(rc,2)=0;
        end
        clear pred_patch_name
    end % if
end

% Now calculate ROC. and then set the threshold to 0.5 and find other
% performance results
test_target=res_table(:,2);
pred=res_table(:,1);

[~,~,~,auc] = perfcurve(test_target,pred,1)

thr=0.5; %threshold for measuring other metrics

predc=pred>thr;
% MCC prec recall fmeas acc
tp=0; fp=0; tn=0; fn=0;
%Measure tp fp tn fn
for i1=1:length(test_target)
    if (test_target(i1)==1)&&(predc(i1)==1), tp=tp+1; end %true positive
    if (test_target(i1)==0)&&(predc(i1)==1), fp=fp+1; end %false positive
    if (test_target(i1)==0)&&(predc(i1)==0), tn=tn+1; end %true negative
    if (test_target(i1)==1)&&(predc(i1)==0), fn=fn+1; end %false negative
end

% Calculate The Matthews correlation coefficient (MCC) or phi coefficient 
sen=tp/(tp+fn)
spe=tn/(tn+fp)
prec=tp/(tp+fp)
acc=(tp+tn)/(tp+tn+fn+fp)
fmeas=2*tp/(2*tp+fp+fn)
mcc=(tp*tn-fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
auc

[sen spe prec acc fmeas mcc auc]

% Command window output on 9 May 2022
% mcc =
%     0.5234
% prec =
%     0.7500
% rec =
%     0.5806
% fmeas =
%     0.6545
% acc =
%     0.7979
% auc =
%     0.8868


