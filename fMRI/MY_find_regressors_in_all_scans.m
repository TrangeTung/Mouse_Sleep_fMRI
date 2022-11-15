function Multi_Regessor = MY_find_regressors_in_all_scans(animal_root_path,path_sub,func_ids,frames,Reg_choices,varargin)

% animal_root_path : 'A:\..\...._1_1\'
% path_sub : Results
% func_ids : 6/7/8
% varargin : {1} experimenttype;
%            {2}{1} raw_rp_file_name;
%            {2}{2} raw_epi_raw_4D_file_name; relative URL
%            {2}{3} raw_epi_pre_4D_file_name; relative URL
%            {3}{1} raw_epi_4D_mask_name; relative URL; Mus_mask.nii//muscle outside the brain
%            {3}{2} raw_epi_4D_mask_name; relative URL; WM_mask.nii
%            {3}{3} raw_epi_4D_mask_name; relative URL; CSF_mask.nii
%            {3}{4} raw_epi_4D_mask_name; relative URL; GS_mask.nii

Multi_Reg = [];%
for ids = 1:numel(Reg_choices)
    str_ram = Reg_choices{ids};
    Reg_choices_str = str_ram(~isstrprop(str_ram,'digit'));
    URL = [animal_root_path,'\',path_sub];
    switch Reg_choices_str
        case {'no'}
            Reg = '';
        case {'rp','rp"','rp2'}
            Reg = MY_get_rp_regessor(varargin{1}{2}{1},Reg_choices{ids},[URL,'\',num2str(func_ids{:}),'\']);
        case {'PCs'}
            Reg = MY_get_img_regressor(varargin{1}{2}{2},varargin{1}{3}{1},Reg_choices{ids},URL);
        case {'WM'}
            Reg = MY_get_img_regressor(varargin{1}{2}{3},varargin{1}{3}{2},Reg_choices{ids},URL);
        case {'CSF'}
            Reg = MY_get_img_regressor(varargin{1}{2}{3},varargin{1}{3}{3},Reg_choices{ids},URL);
        case {'GS'}
            Reg = MY_get_img_regressor(varargin{1}{2}{3},varargin{1}{3}{4},Reg_choices{ids},[URL,'\',num2str(func_ids{:}),'\']);
    end
    Multi_Reg = [Multi_Reg Reg];
end
Multi_Reg(1:min(frames(:))-1,:) = [];
cd([animal_root_path,'\',path_sub,'\',num2str(func_ids{:})])
for k = 1:size(Multi_Reg,2)
   RAM =  Multi_Reg(:,k);
   Multi_Reg(:,k) = (RAM-mean(RAM))/std(RAM); 
end
save 'Multi_Regessor.txt' 'Multi_Reg' '-ascii'
switch varargin{1}{1}
    case 1
        Multi_Regessor = [animal_root_path,'\',path_sub,'\',num2str(func_ids{:}) '\Multi_Regessor.txt'];
    case 2
        Multi_Regessor = [ones(size(Multi_Reg,1),1),Multi_Reg];
end

end

function Reg = MY_get_rp_regessor(rp_text,para,URL)

rp_raw = load([URL,rp_text]);
switch para
    case {'rp'}
        Reg = rp_raw;
    case {'rp"'}
        Reg = [rp_raw,rp_raw - [zeros(1,6);rp_raw(1:end-1,:)]];
    case {'rp2'}
        rp_dev = rp_raw - [zeros(1,6);rp_raw(1:end-1,:)];
        Reg = [rp_raw rp_dev].^2;
end
end

function Reg = MY_get_img_regressor(img_file,mask_file,para,URL)

fMRI_data = spm_read_vols(spm_vol([URL,'\',img_file]));
%lmask = spm_read_vols(spm_vol([URL,'\..\',mask_file]));
lmask = spm_read_vols(spm_vol([URL,'\..\',mask_file]));
para_str = (para(~isstrprop(para,'digit')));

switch para_str
    case {'WM','CSF'}
        PC_num = 1;
        lmask = MY_erode_mask(lmask,3);
        fMRI_RAM = fmask(fMRI_data,lmask);
        Reg = MY_get_principle_component(fMRI_RAM,PC_num);
    case {'GS'}
        fMRI_RAM = fmask(fMRI_data,lmask);
        Reg = mean(fMRI_RAM,1);
        Reg = Reg(:);
    case {'PCs'} 
        PC_num = str2double(para(isstrprop(para,'digit')));
        lmask = MY_erode_mask(lmask,3);
        fMRI_RAM = fmask(fMRI_data,lmask);
        Reg = MY_get_principle_component(fMRI_RAM,PC_num);
end
end


function PCs = MY_get_principle_component(fMRI_RAM,PC_num)

fMRI_noise_ready_PCA = (fMRI_RAM-mean(fMRI_RAM,2))./std(fMRI_RAM,0,2);
[~,score,latent,~]= pca(fMRI_noise_ready_PCA','algorithm','svd');
latent_test = sort(latent,'descend');
if sum(latent_test-latent) ~= 0
    error('Something was wrong during PCA');
else
    PCs = score(:,1:PC_num);
end
end

function lmask = MY_erode_mask(mask,radius_erode)

SE = strel('square',radius_erode);
mask_RAM = reshape(mask,[size(mask,1) size(mask,2)*size(mask,3)]);
mask_RAM =  imerode(mask_RAM,SE);
lmask = reshape(mask_RAM,size(mask));
end
