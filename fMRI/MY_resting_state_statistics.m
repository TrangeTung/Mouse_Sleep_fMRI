function MY_resting_state_statistics(animal_root_path,path_sub,func_ids,frames,Reg_choices,defined_rest,result_1st)

group_scan = defined_rest.Nscans;
switch group_scan
    case 'individual'
        individual_resting_correlation(animal_root_path,path_sub,func_ids,frames,Reg_choices,defined_rest,result_1st);
    case 'Allscans'
        allscans_resting_correlation(animal_root_path,path_sub,func_ids,frames,Reg_choices,defined_rest,result_1st);
end

%% Fig2PPTX
seedname = defined_rest.seedname;
EPI_folder = func_ids{:};
dest = [animal_root_path,'\Functions\rsfMRI'];
fig2PPTX(dest,EPI_folder,seedname,group_scan);


end

function individual_resting_correlation(animal_root_path,path_sub,func_ids,frames,Reg_choices,defined_rest,result_1st)
%function individual_resting_correlation(animal_root_path,func_ids,frames,Reg_choices,defined_rest,result_1st)

seedname   = defined_rest.seedname;
seedpath   = defined_rest.seedpath;
Tem_name   = defined_rest.Tem_name;
filename   = defined_rest.filename;
filepath   = defined_rest.filepath;
Gmask_name = defined_rest.denoi_profile;
slice      = result_1st.slice;
colorbar   = result_1st.colorbar;

EPI_folder = func_ids{:};
dest = [animal_root_path,'\Functions\rsfMRI'];
lmask = spm_read_vols(spm_vol(Gmask_name));
co_ef = zeros([size(lmask),size(seedname,1),size(seedname,2),length(EPI_folder)]);
for k = 1:length(EPI_folder)
    dest_cc = [dest,'\',num2str(EPI_folder(k))];
    if ~exist(dest_cc,'dir');mkdir(dest_cc); end
    cd(dest_cc);
    if exist(['r',filename,'.nii'],'file') 
        fhead = spm_vol(['r',filename,'.nii,1']);ram = char(fhead.descrip); 
    else
        ram = char();
    end
    try
        Reg_ram = transpose(char(Reg_choices));
        Center_heart_test = prod(strtrim(ram) == strtrim(Reg_ram(:)'));
        if Center_heart_test ==1;fprintf('4D images have existed.');end
        all_img = spm_read_vols(spm_vol(['r',filename,'.nii']));
    catch
        cd ([filepath '\' num2str(EPI_folder(k))]);
        Segments = MY_search_bruker_method('Segments',EPI_folder(k),animal_root_path);
        EPI_TR = MY_search_bruker_method('EPI_TR',EPI_folder(k),animal_root_path)/1000*Segments;
        %% regressors
        varargin = {2;{'rp_s2dseq.txt';'2dseq.nii';[filename,'.nii']};...
            {'Mus_mask.nii';'WM_mask.nii';'CSF_mask.nii';'GS_mask.nii'}};
        RegBasFuc = MY_find_regressors_in_all_scans(animal_root_path,path_sub,{EPI_folder(k)},frames,Reg_choices,varargin);
        
        head = spm_vol([filename,'.nii']);
        all_img = spm_read_vols(head);
        data_ready_regress = fmask(all_img,lmask);
        data_ready_regress(isnan(data_ready_regress)) = 0;
        for iii = 1:size(data_ready_regress,1)
            [Beta,~,Residual] = regress(squeeze(data_ready_regress(iii,:))',RegBasFuc);
            data_ready_regress(iii,:) = Residual + Beta(1);
        end
        data_regressout = funmask(data_ready_regress,lmask);
        %% bandpass filter
        all_img = Smooth_temporal_f(data_regressout,lmask, EPI_TR);
        ram = transpose(char(Reg_choices));
        for n = 1:numel(head)
           head(n).fname = ['r',filename,'.nii'];
           head(n).descrip = ram(:)';%
        end
        cd(dest_cc);
        spm_write_vol_4D(head,all_img);
    end
    %% Normalize rs data
    all_img = fmask(all_img,lmask);
    all_img = all_img - repmat(mean(all_img,2),[1,size(all_img,2)]);
    all_img = funmask(all_img,lmask);
    
    for i = 1:size(seedname,1)
        for j = 1:size(seedname,2)
            fMRI_data = all_img;
            seedmask = spm_read_vols(spm_vol([seedpath,'\',seedname{i,j},'.nii']));
            seedseque = mean(fmask(fMRI_data,seedmask),1);
            
            fMRI_data(:,:,:,1:min(frames(:))-1) = [];
            seedseque(1:min(frames(:))-1) = [];
            TC(i,j,k,:) = seedseque;
            
            f_signal = fmask(fMRI_data,lmask);
            coef_f_seed = corr(f_signal',seedseque');
            co_ef(:,:,:,i,j,k) = funmask(coef_f_seed,lmask);
            head = spm_vol(Tem_name);
            mapname = ['CCmap_' seedname{i,j} '_' num2str(EPI_folder(k))];
            head.fname = ['CCmap_' seedname{i,j} '_' num2str(EPI_folder(k)) '.nii'];
            cd(dest_cc);
            head.pinfo(1:2)=[0;1];
            spm_write_vol(head,flip(squeeze(co_ef(:,:,:,i,j,k)),1));
            
            Colormap(   'statfile',[dest_cc,'\',head.fname],...
                'bgmfile',Tem_name,...
                'denoi_profile',Gmask_name,...
                'slice',slice,...
                'bar_value',colorbar,...
                'dest',dest_cc,...
                'mapname',mapname,...
                'cluster',10);
            
            epi_Tem_name = [filepath '\' num2str(EPI_folder(k)) '\' filename '.nii,1'];
            Colormap(   'statfile',[dest_cc,'\',head.fname],...
                'bgmfile',epi_Tem_name,...
                'denoi_profile',Gmask_name,...
                'slice',slice,...
                'bar_value',colorbar,...
                'dest',dest_cc,...
                'mapname',[mapname,'_epi'],...
                'cluster',10  );
        delete([dest_cc,'\',mapname,'_epi_nothres.tif'])
        end
    end
end

%% XLSwrite
nPointFFT=1024;
Segments = MY_search_bruker_method('Segments',EPI_folder(k),animal_root_path);
EPI_TR = MY_search_bruker_method('EPI_TR',EPI_folder(k),animal_root_path)/1000*Segments;
fs = 1/EPI_TR;
f = 0:fs/2/nPointFFT:fs/2;
count = 0; CC_excel = 'CC_Value.xlsx';
for i = 1:size(seedname,1)
    outexcel = strcat(seedname{i,1}(1:end-1),'.xlsx');
    for j = 1:size(seedname,2)
        count = count+1;
        for k = 1:numel(EPI_folder)
            cd(dest);
            %% Long-range connectivity
            seed_mask = spm_read_vols(spm_vol([seedpath,'\',seedname{i,rem(j,2)+1},'.nii']));
            seed_RAM = reshape(seed_mask,[size(seed_mask,1) size(seed_mask,2)*size(seed_mask,3)]);
            SE = strel('square',4);
            seed_RAM = imdilate(seed_RAM,SE);
            seed_mask = reshape(seed_RAM,size(seed_mask));
            
            CC_map = spm_read_vols(spm_vol([dest,'\',num2str(EPI_folder(k)),'\','CCmap_' seedname{i,j} '_' num2str(EPI_folder(k)) '.nii']));
            CC_seed = CC_map.*seed_mask;
            CC_Value = mean(CC_seed(CC_seed(:)~=0&~isnan(CC_seed(:))));
            
            xlswrite(outexcel,EPI_folder(k),'CC_Value',strcat(char(double('A')+k),'1'));
            xlswrite(outexcel,CC_Value,'CC_Value',strcat(char(double('A')+k),num2str(1+j)));
            xlswrite(CC_excel,EPI_folder(k),'CC_Value',strcat(char(double('A')+k),'1'));
            xlswrite(CC_excel,{seedname{i,j}},'CC_Value',strcat(char(double('A')),num2str(1+count)));
            xlswrite(CC_excel,CC_Value,'CC_Value',strcat(char(double('A')+k),num2str(1+count)));
            % ------------------------------------------------
            % delete it
            %% Local connectivity
            Lseed_mask = spm_read_vols(spm_vol([seedpath,'\',seedname{i,j},'.nii']));
            Lseed_RAM = reshape(Lseed_mask,[size(Lseed_mask,1) size(Lseed_mask,2)*size(Lseed_mask,3)]);
            SE = strel('disk',4);
            Lseed_RAM = imdilate(Lseed_RAM,SE);
            Lseed_mask = reshape(Lseed_RAM,size(Lseed_mask));
            
            CC_map = spm_read_vols(spm_vol([dest,'\',num2str(EPI_folder(k)),'\','CCmap_' seedname{i,j} '_' num2str(EPI_folder(k)) '.nii']));
            LCC_seed = CC_map.*Lseed_mask;
            Local_CC_Value = mean(LCC_seed(LCC_seed(:)~=0&~isnan(LCC_seed(:))));
            
            xlswrite(outexcel,EPI_folder(k),'Local_CC_Value',strcat(char(double('A')+k),'1'));
            xlswrite(outexcel,Local_CC_Value,'Local_CC_Value',strcat(char(double('A')+k),num2str(1+j)));
            xlswrite(CC_excel,EPI_folder(k),'Local_CC_Value',strcat(char(double('A')+k),'1'));
            xlswrite(CC_excel,{seedname{i,j}},'Local_CC_Value',strcat(char(double('A')),num2str(1+count)));
            xlswrite(CC_excel,Local_CC_Value,'Local_CC_Value',strcat(char(double('A')+k),num2str(1+count)));
            % -------------------------------------------------
            %% Time course
            TC_seed = squeeze(TC(i,j,k,:));
            Sheet = strcat('TC_',seedname{i,j});
            Range = strcat(char(double('A')+k),'2');
            xlswrite(outexcel,TC_seed,Sheet,Range);
            xlswrite(outexcel,EPI_folder(k),Sheet,strcat(char(double('A')+k),'1'));
            %% pwelch. calculate the power spectrum
            TC_seedN = (TC_seed - mean(TC_seed(:)))/std(TC_seed(:));
            [pxx,fpwelch] = pwelch(TC_seedN,[],[],f,fs);
            Sheet = strcat('pwelch_',seedname{i,j});
            Range = strcat(char(double('A')+k),'2');
            xlswrite(outexcel,pxx',Sheet,Range);
            xlswrite(outexcel,EPI_folder(k),Sheet,strcat(char(double('A')+k),'1'));
            %% fft
            [seedfft,ffft] = periodogram(TC_seedN,[],f,fs);
            Sheet = strcat('fft_',seedname{i,j});
            Range = strcat(char(double('A')+k),'2');
            xlswrite(outexcel,seedfft',Sheet,Range);
            xlswrite(outexcel,EPI_folder(k),Sheet,strcat(char(double('A')+k),'1'));
        end
        xlswrite(outexcel,{seedname{i,j}(1:end)},'CC_Value',strcat('A',num2str(1+j)));
        xlswrite(outexcel,(1:numel(TC_seed))',strcat('TC_',seedname{i,j}),'A2');
        xlswrite(outexcel,fpwelch',strcat('pwelch_',seedname{i,j}),'A2');
        xlswrite(outexcel,ffft',strcat('fft_',seedname{i,j}),'A2');
    end
end

end

function allscans_resting_correlation(animal_root_path,path_sub,func_ids,frames,Reg_choices,defined_rest,result_1st)

seedname   = defined_rest.seedname;
seedpath   = defined_rest.seedpath;
Tem_name   = defined_rest.Tem_name;
filename   = defined_rest.filename;
filepath   = defined_rest.filepath;
Gmask_name = defined_rest.denoi_profile;
EPI_folder = func_ids{:};

slice      = result_1st.slice;
colorbar   = result_1st.colorbar;

for i = 1:size(seedname,1)
    for j = 1:size(seedname,2)
        dest = [animal_root_path,'\Functions\rsfMRI'];
        try
            cd([dest,'\',num2str(EPI_folder(1))])
        catch
            individual_resting_correlation(animal_root_path,path_sub,func_ids,frames,Reg_choices,defined_rest,result_1st)
        end
        head = spm_vol(Tem_name);
        template = spm_read_vols(head);
        lmask = (template>0.01);
        seedmask = spm_read_vols(spm_vol([seedpath,'\',seedname{i,j},'.nii']));
        co_ef = zeros([size(seedmask) length(EPI_folder)]);
        
        for k = 1:length(EPI_folder)
            cd([dest,'\',num2str(EPI_folder(k))])
            co_ef(:,:,:,k) =spm_read_vols(spm_vol(['CCmap_' seedname{i,j} '_' num2str(EPI_folder(k)) '.nii']));
        end
        %% mean CCMap for all trails
        dest_cc = [animal_root_path,'\Functions\rsfMRI\Mean'];
        if ~exist(dest_cc,'dir');mkdir(dest_cc);end
        mean_co_ef = mean(co_ef,4);
        mapname = ['CCmap_' seedname{i,j} '_mean'];
        head.fname = ['CCmap_' seedname{i,j} '_mean.nii'];
        cd(dest_cc);
        spm_write_vol(head,squeeze(mean_co_ef));
        
        Colormap(   'statfile',[dest_cc,'\',head.fname],...
                    'bgmfile',Tem_name,...
                    'denoi_profile',Gmask_name,...
                    'slice',slice,...
                    'bar_value',colorbar,...
                    'dest',dest_cc,...
                    'mapname',mapname,...
                    'cluster',10  )
        
        epi_Tem_name = [filepath '\' num2str(EPI_folder(end)) '\' filename '.nii,1'];
        Colormap(   'statfile',[dest_cc,'\',head.fname],...
                    'bgmfile',epi_Tem_name,...
                    'denoi_profile',Gmask_name,...
                    'slice',slice,...
                    'bar_value',colorbar,...
                    'dest',dest_cc,...
                    'mapname',[mapname,'_epi'],...
                    'cluster',10  )
        delete([dest_cc,'\',mapname,'_epi_nothres.tif'])
    end
end

end

function fig2PPTX(path,folder,seedname,Nscans)
cd([path '\' num2str(folder(1))]);
tif_dir = dir('*.tif');
tif_img = imread(tif_dir(1).name);
PPT_length = 16;
PPT_width = 9;
gap = 0.05;
ppt_fig_ud = PPT_length/20*9*size(tif_img,1)/size(tif_img,2)+gap;
ppt_fig_lr = PPT_length/20*9;
cd(path);cd ..;
try
    exportToPPTX('close');
    exportToPPTX('open', 'rsfMRI.pptx');
catch
    exportToPPTX('new','Dimensions',[PPT_length PPT_width]);
end

switch Nscans
    case 'Allscans'
        %% Mean image
        cd([path '\Mean']);
        exportToPPTX('addslide');
        exportToPPTX('addtext', 'Mean C.C.Map of all trials', 'Position', [PPT_length/4 0 PPT_length/2 0.5],'HorizontalAlignment','center','FontWeight','bold');
        for j = 1:size(seedname,1)
            %% left
            exportToPPTX('addtext',seedname{j,1},'Position',[0 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
            img = imread(strcat('CCmap_',seedname{j,1},'_mean.tif'));
            exportToPPTX('addpicture',img,'Position',[PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
            %% right
            exportToPPTX('addtext',seedname{j,2},'Position',[PPT_length/2 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
            img = imread(strcat('CCmap_',seedname{j,2},'_mean.tif'));
            exportToPPTX('addpicture',img,'Position',[PPT_length/2+PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
        end
        
        %% Mean image ¡ª¡ªepi template
        cd([path '\Mean']);
        exportToPPTX('addslide');
        exportToPPTX('addtext', 'Mean C.C.Map of all trials', 'Position', [PPT_length/4 0 PPT_length/2 0.5],'HorizontalAlignment','center','FontWeight','bold');
        for j = 1:size(seedname,1)
            %% left
            exportToPPTX('addtext',seedname{j,1},'Position',[0 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
            img = imread(strcat('CCmap_',seedname{j,1},'_mean_epi.tif'));
            exportToPPTX('addpicture',img,'Position',[PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
            %% right
            exportToPPTX('addtext',seedname{j,2},'Position',[PPT_length/2 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
            img = imread(strcat('CCmap_',seedname{j,2},'_mean_epi.tif'));
            exportToPPTX('addpicture',img,'Position',[PPT_length/2+PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
        end
    case 'individual'
        %% CC for different trials
        cd(path);cd ..;
        %% single trial
        for i = folder(:)'
            cd([path '\' num2str(i)]);
            exportToPPTX('addslide');
            exportToPPTX('addtext', strcat('C.C.Map of trial_',num2str(i)), 'Position', [PPT_length/4 0 PPT_length/2 0.5],'HorizontalAlignment','center','FontWeight','bold');
            for j = 1:size(seedname,1)
                %% left
                exportToPPTX('addtext',seedname{j,1},'Position',[0 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
                img = imread(strcat('CCmap_',seedname{j,1},'_',num2str(i),'.tif'));
                exportToPPTX('addpicture',img,'Position',[PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
                %% right
                exportToPPTX('addtext',seedname{j,2},'Position',[PPT_length/2 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
                img = imread(strcat('CCmap_',seedname{j,2},'_',num2str(i),'.tif'));
                exportToPPTX('addpicture',img,'Position',[PPT_length/2+PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
            end
        end
        
        %% CC for different trials __ epi template
        cd(path);cd ..;
        %% single trial
        for i = folder(:)'
            cd([path '\' num2str(i)]);
            exportToPPTX('addslide');
            exportToPPTX('addtext', strcat('C.C.Map of trial_',num2str(i)), 'Position', [PPT_length/4 0 PPT_length/2 0.5],'HorizontalAlignment','center','FontWeight','bold');
            for j = 1:size(seedname,1)
                %% left
                exportToPPTX('addtext',seedname{j,1},'Position',[0 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
                img = imread(strcat('CCmap_',seedname{j,1},'_',num2str(i),'_epi.tif'));
                exportToPPTX('addpicture',img,'Position',[PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
                %% right
                exportToPPTX('addtext',seedname{j,2},'Position',[PPT_length/2 0.5+ppt_fig_ud*(j-1) PPT_length/20 1],'Scale','noscale','FontSize',10,'HorizontalAlignment','right','FontWeight','bold');
                img = imread(strcat('CCmap_',seedname{j,2},'_',num2str(i),'_epi.tif'));
                exportToPPTX('addpicture',img,'Position',[PPT_length/2+PPT_length/20 0.5+ppt_fig_ud*(j-1) ppt_fig_lr ppt_fig_ud-gap],'Scale','noscale');
            end
        end
        
end
cd(path);cd ..;
exportToPPTX('save', 'rsfMRI.pptx');
exportToPPTX('close');
end
