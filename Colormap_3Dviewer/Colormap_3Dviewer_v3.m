function [f1_rgb]=Colormap_3Dviewer_v3(Func3D,Temp3D,bar,codepath,frame,slice_perc)

% version 3 

[pathstr,name,ext]=fileparts(Func3D);
func = fullfile(pathstr,['n_',name,ext]);
if ~exist(func,'file')
    ref{1,1} = fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20_I2.nii');
    source{1,1} = Temp3D;
    all_func = MY_select_file_for_SPM(pathstr,['^',name,ext],[1 Inf]);
    OldNormalize_mlb = MY_get_default_oldnormalize_batch_struct(ref, source, all_func);
    OldNormalize_mlb{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'n_';
    OldNormalize_mlb{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
    F = spm_figure('GetWin');
    spm_jobman('run',OldNormalize_mlb);
    hgexport(figure(F), fullfile(pwd, strcat('oldnormalize')), hgexport('factorystyle'), 'Format', 'tiff');
    
    all_func = MY_select_file_for_SPM(pathstr,['^n_',name,ext],[1 Inf]);
    Smooth_mlb = MY_get_default_smooth_batch_struct(all_func,8);
    disp('Start to process Smooth!');
    spm_jobman('run',Smooth_mlb);

end

tname = fullfile(codepath,'Colormap_3Dviewer','average_template_50.nii');
L = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Label_Mouse_X20.nii')));
timg = spm_read_vols(spm_vol(tname));
timg = flip(permute(timg,[3 2 1]),3);
lmask = timg>0;

fimg_ = spm_read_vols(spm_vol([func,',',num2str(frame)])); 
fimg_ = fimg_/2+flip(fimg_,1)/1.5;
fimg_(fimg_>bar(end))=bar(end);
fimg_(fimg_<bar(1))=bar(1);
fimg_(isnan(fimg_))=0;
fimg = imresize3(fimg_,size(timg));
fimg = medfilt3(fimg,[1 1 1]*9);
% fimg = smooth3(fimg,'box',11);
fimg(~lmask)=nan;
fimg(isnan(fimg))=0;

timg(timg>200)=200;
BkVol = timg; 
%BkVol = imresize3(BkVol,3);
BkVol = smooth3(BkVol,'box',3);

LbVol = fimg;
LbVol(:,1:2,:) = [];
LbVol = imresize3(LbVol,size(fimg));
%LbVol = imresize3(LbVol,3);
LbVol = smooth3(LbVol,'box',5);

Nslice = size(BkVol,3);
Cut = round(Nslice-Nslice*slice_perc);
BkVol(:,:,end-Cut:end)=0;
LbVol(:,:,end-Cut:end)=0;

Ls=LbVol*0; Mis=zeros(121,1);
for loop=1:61
    if loop~=61
        Dif1=(bar(2)-bar(1)) / 60;
        Dif2=(bar(4)-bar(3)) / 60;

        Cs1 = LbVol>bar(1)+(loop-1)*Dif1 & LbVol<=bar(1)+(loop)*Dif1;
        Cs2 = LbVol<bar(4)-(loop-1)*Dif2 & LbVol>=bar(4)-(loop)*Dif2;
        
        if isempty(find(Cs1(:)==1));Mis(loop)=1;end
        if isempty(find(Cs2(:)==1));Mis(121-loop)=1;end

        Ls(Cs1) = loop;
        Ls(Cs2) = 121-loop;
    end
    if loop==61
        Cs = LbVol>bar(2) & LbVol<bar(3);
        Ls(Cs) = 0;
    end
end
Ls(LbVol<=bar(1))=1;
Ls(LbVol>=bar(4))=120;
Ls(isnan(Ls))=0;


MyMap = [ones(1,3);flipud(winter(60));autumn(60)];
% MyMap = [ones(1,3);lines(120)];
domap = [ones(128,1),(63.5:0.5:127)'/127,(0:127)'/127];
gdmap = [40/255+(255-40)/255*(0:127)'/127,125/255+(255-125)/255*(0:127)'/127,184/255+(255-184)/255*(0:127)'/127];
% MyMap = [ones(1,3);gdmap(1:60,:);flipud(domap(1:60,:));];


MyMap(Mis==1,:)=[];
BiVis = [false;true(120,1)];
BiVis(Mis==1,:)=[];BiVis(1)=false;
LbOpcy = ones(numel(BiVis),1)*0.9;


config.CameraPosition=[-2 2.5 3];
config.CameraUpVector = [-1 0 0];
config.CameraTarget = [0 0 0];
config.CameraViewAngle = 15;
config.BackgroundColor = [0 0 0]+1;
config.ShowIntensityVolume = 1;
config.LabelColor = MyMap;
config.LabelVisibility = BiVis;
config.LabelOpacity = LbOpcy;
config.VolumeOpacity = 0.5;
config.VolumeThreshold = 0.5;

F1=figure;
config.CameraPosition=[-2,2,2];
config.CameraUpVector = [-1 0 0];
labelvolshow(Ls,BkVol,config);
f1_rgb = getframe(F1);


end



