%% preparation
% According to your experimental parameters and certain file direction, Change these bellow.
clear all;clc
codepath = 'C:\Users\yuyalin\Desktop\rsfMRI';
cd(codepath)
addpath(genpath(codepath));
%% load
Animal_path{1}='';  %fmri data
epath{1}='';        % state label

iDataNum = size(Animal_path,2);
cd(Animal_path{1});
template = 'D:\lianglab\code\fmri\mouse\Template_Mouse_v38.nii';
Animaltype = 'mouse';
%% SPM para
% SPM_para = struct()
sliceorder = 'IA';
% sliceorder = input('Please Input Slice order, eg: IA ','s' );
Reg_choices = 'rp';
% Reg_choices = input('Please Input regression method, eg: rp ', 's');
PEdirection = 2;
% PEdirection = str2num(input('Please Input PE direction, eg: 2 ','s' ));
%% other para
pvalue = 0.01;
smooth_para = [8 8 8];
%% regressor choices
%Reg_choices = {'rp';'rp"';}; % For task-based analysis
%'GS';
% no     --->  no regression
% rp    --->  6 head motion parameters regression [rp]
% rp"   --->  the deviation of 6 head motion parameters regression [rp']
% rp2   --->  the square of rp and rp" [rp.^2 rp'.^2]
% XPCs  --->  X Principal Components outside the brain
% WM    --->  the Principal Components in white matter mask "WM_mask.nii";   // storage in the folder '..\Results\'
% CSF   --->  the Principal Components in CSF mask "CSF_mask.nii";           // storage in the folder '..\Results\'
% GS    --->  the global signal in global mask "GS_mask.nii";                // storage in the folder '..\Results\'
% the PC was extraction automatically by c

inEPI_number = 0;
spm('defaults', 'FMRI');
set(spm('CreateIntWin','off'),'Visible','on');
for number = [1]%1:iDataNum
    datapath = [Animal_path{number},'\'];
    % Read every scan method
    subdir =  dir( datapath );
    EPIrecord = 0;
    EPI_folder = Animal_EPI_folder{number};
    for scan = 3:length(subdir)
        if subdir(scan).isdir
            judge_folder = isstrprop(subdir(scan).name,'digit');
            if judge_folder(1) == 1 % is scan data 
                % read method to confirm seq type
                methodpath = [subdir(scan).folder,'\',subdir(scan).name,'\method'];
                seq_name = read_method(methodpath);
                switch seq_name
                    case 'RARE'
                    RARE = subdir(scan).name;  
                    case 'FieldMap'
                    Fieldmap  = subdir(scan).name;  
%                     case 'EPI'
%                     EPI_folder = str2num(subdir(scan).name);
                end  
            end
        end
    end
    Segments = MY_search_bruker_method('Segments',EPI_folder(1),datapath);
    EPI_TR = MY_search_bruker_method('EPI_TR',EPI_folder(1),datapath)/1000*Segments;

    for flag_stage = [1:7]
            
        if flag_stage == 1
          %% Slicetiming
                Nslice = MY_search_bruker_method('Nslice',EPI_folder(1),datapath);
                all_func = MY_find_images_in_all_scans(datapath,'Results',{EPI_folder},'^m2dseq','.nii',[1 Inf],'separate_cells');
                slicetiming_mlb = MY_get_default_slicetiming_batch_struct(all_func,Nslice,EPI_TR,sliceorder);
                disp('Start to process Slicetiming!')
                spm_jobman('run',slicetiming_mlb);
        end
        if flag_stage == 2
          %% Realignment
            all_func = MY_find_images_in_all_scans(datapath,'Results',{EPI_folder},'^sm2dseq','.nii',[1 Inf],'separate_cells');
            realign_mlb = MY_get_default_realign_batch_struct(all_func);
            F = spm_figure('GetWin');
            disp('Start to process realignment !')
            spm_jobman('run',realign_mlb);
            hgexport(figure(F), fullfile([datapath,'Results\'],strcat('realign')), hgexport('factorystyle'), 'Format', 'tiff');
            clear realign_mlb all_func;      
        end   
    
        if flag_stage == 3
          %% Func2T2 Coregistration
            ref{1,1} = [datapath 'Results\' num2str(EPI_folder(1)) '\rrsm2dseq.nii,1'];
            source{1,1} = [datapath 'Results\T2\T2_m.nii,1'];
            Func2T2W_mlb = MY_get_default_coreg_batch_struct(ref, source, {''});
            disp('Start to process Func2T2W coregistration!');
            F = spm_figure('GetWin');
            spm_jobman('run',Func2T2W_mlb);
            hgexport(figure(F), fullfile([datapath 'Results\'], 'coreg'), hgexport('factorystyle'), 'Format', 'tiff');

            clear Func2T2W_mlb other ref source;          
        end
        
        if flag_stage == 4
          %% T22Template coregistration
            ref{1,1} = template;
            source{1,1} = [datapath 'Results\T2\cT2_m.nii,1'];
            all_func = MY_find_images_in_all_scans(datapath,'Results',{EPI_folder(:)},'^rrsm2dseq','.nii',[1 Inf],'all_mixed');
            all_func = [all_func;[datapath 'Results\T2\cT2_m.nii,1']];
            OldNormalize_mlb = MY_get_default_oldnormalize_batch_struct(ref, source, all_func);
            disp('Start to process OldNormalize!');
            F = spm_figure('GetWin');
            spm_jobman('run',OldNormalize_mlb);
            hgexport(figure(F), fullfile([datapath 'Results\'], strcat('oldnormalize')), hgexport('factorystyle'), 'Format', 'tiff');
            clear OldNormalize_mlb other ref source;
        end
        
        if flag_stage == 5
          %% smooth_space
            all_func = MY_find_images_in_all_scans(datapath,'Results',{EPI_folder},'^Lnrrsm2dseq','.nii',[1 Inf],'all_mixed');
            Smooth_mlb = MY_get_default_smooth_batch_struct(all_func,smooth_para);
            disp('Start to process Smooth!');
            spm_jobman('run',Smooth_mlb);
            clear Smooth_mlb;
        end
                
        if flag_stage == 6
          %% non-brain PCA
                hh=[datapath  'Results\' num2str(EPI_folder(1))];
                
                fMRI_data = spm_read_vols(spm_vol([datapath,'Results\' num2str(EPI_folder(1)) '\2dseq.nii']));
                lmask = spm_read_vols(spm_vol([datapath,'Rsults\' num2str(EPI_folder(1)) '\EPI_mask.nii']));
                cd('D:\software\matlab\toolbox\stats\stats');
                %PC = MY_get_principle_component_out_of_brian(fMRI_data,lmask,80);
                coeff=MY_get_principle_component_out_of_brian(fMRI_data,lmask,80);
                cd (hh)
                save 'coeff.txt' -ascii coeff;
                %PC1=PC(:,1:40);
%                 load rp_sm2dseq.txt
%                 rp_sm2dseq(1,7:12) = 0;
%                 rp_sm2dseq(2:size(rp_sm2dseq,1),7:12) = rp_sm2dseq(2:end,1:6) - rp_sm2dseq(1:end-1,1:6);
                
%                 save 'PC_all.txt' -ascii PC;
%                 save 'PC40.txt' -ascii PC1;

 %               save 'rp.txt' -ascii rp_sm2dseq; 
        end       
           
        if flag_stage == 7            
                cd(epath{number});
                load('aw.mat');
                duration=cell(1,3);
                onset=cell(1,3);
                duration{1}=wake_dura;onset{1}=onsets;
                duration{3}=nrem_dura;onset{3}=nrem_onsets;
                duration{4}=rem_dura;onset{4}=ronsets;
                
                ww=cell(1,3);
                ww{1}=1;
                ww{2}=[0 1];
                ww{3}=[0 0 1];
                
                folder = EPI_folder;
              %%                
                Reg_choices = {'rp';};
                result_1st = struct('weights',ww,'slice',1:38,'template',template,'FDR_pvalue',pvalue,'colorbar',colorbar);                
                defined_1st = struct('Nscans','Allscans','filename','sLnrrsm2dseq','duration',duration,'onset',onset);
              %% change spmhrf.m
%                 if Animaltype == 'mouse'
%                     cd(SPM_path)
%                     eval(['!rename' ,   '/,spm_hrf.m', ',spm_hrf_back.m']);
%                     cd(codepath)
%                     copyfile 'spm_hrf.m' SPM_path;
%                 end

                MY_task_state_statistics(datapath,'Results',{folder(:)},[1 Inf],Reg_choices,defined_1st,result_1st);
                
                % ---------- rename the animal-wise folder ----------
%                 cd([datapath '\Functions\tsfMRI']);
%                 if exist('Allscans','dir');rmdir('Allscans','s');end
%                 mkdir Allscans
                % ---------- delete it after this protocol ----------
        end
         if flag_stage == 8
             
                folder = EPI_folder;
                cd([datapath  'Results\' num2str(folder) ]);
                load rp_sm2dseq.txt
                load PC40.txt
                 rp_sm2dseq(1,7:12) = 0;
                 rp_sm2dseq(2:size(rp_sm2dseq,1),7:12) = rp_sm2dseq(2:end,1:6) - rp_sm2dseq(1:end-1,1:6);
                header = spm_vol('rsm2dseq.nii');
                Image_4d = spm_read_vols(header);
                mask = spm_read_vols(spm_vol('EPI_mask.nii'));%D:\lianglab\code\fmri\mouse\mouse_ln\test_Template.nii
                flat_Image_4d = fmask(Image_4d,mask);
           
                
                regress_Image = zeros(size(flat_Image_4d,1),size(flat_Image_4d,2));
                for voxelNum = 1 : size(flat_Image_4d,1)
                    [b,bint,r] = regress(flat_Image_4d(voxelNum,:)',[ones(size(rp_sm2dseq,1),1),rp_sm2dseq,PC40]);%,PC,]
                    regress_Image(voxelNum,:) = r'+b(1);
                end
                
                regress_Image_4d = funmask(regress_Image,mask);
                for i = 1:numel(header)
                    header(i).fname = ['r',header(i).fname];
                end
                spm_write_vol_4D(header,regress_Image_4d);
                clear regress_Image_4d
                
                
                clear regress_Image
                clear Image_4D
                clear flat_Image_4d
         end
    end
end
%%
scans=cell(46,1);
for i=[1:46]
    scans{i}=[Animal_path{i},'\activation_L_1\con_0002.nii'];
end
dest=['C:\Users\yuyalin\Desktop\acti_0307'];
one_sample_ttest_mlb = MY_2nd_level_analysis_flexible_factorial_design_batch_struct(dest,scans,session_num,subj_num);
F = spm_figure('GetWin');
spm_jobman('run',one_sample_ttest_mlb);
hgexport(figure(F), fullfile(dest,strcat('block design')), hgexport('factorystyle'), 'Format', 'tiff');
estimate_mlb = MY_1st_level_analysis_estimate_batch_struct([dest '\SPM.mat']);
spm_jobman('run',estimate_mlb);
results_mlb = MY_1st_level_analysis_results_batch_struct([dest '\SPM.mat'],'tsfMRI',ones(22,1)/22,'none');
spm_jobman('run',results_mlb);
