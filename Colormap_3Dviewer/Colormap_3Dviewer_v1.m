function [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,bar,codepath,frame)

% version 1 : only for no-threshold map view

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
    Smooth_mlb = MY_get_default_smooth_batch_struct(all_func);
    disp('Start to process Smooth!');
    spm_jobman('run',Smooth_mlb);

end

tname = fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii');
L = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Label_Mouse_X20.nii')));
timg = spm_read_vols(spm_vol(tname));
lmask = timg>0;

fimg_ = spm_read_vols(spm_vol([func,',',num2str(frame)]));
fimg_ = fimg_/2+flip(fimg_,1)/2;
fimg_(fimg_>bar(2))=bar(2);
fimg_(fimg_<bar(1))=bar(1);
fimg_(isnan(fimg_))=0;
fimg = imresize3(fimg_,size(timg));
fimg = medfilt3(fimg,[1 1 1]*9);
% fimg = smooth3(fimg,'box',11);
fimg(~lmask)=nan;

Ns=L*0;Ns(size(Ns,1)/2+1:end,:,:)=1;
fimg(Ns&L<=416)=nan;
fimg(L>=1100)=nan;

Vx = fimg(~isnan(fimg));
if min(Vx(:))>bar(1);fimg(10,10,10)=bar(1);end
if max(Vx(:))<bar(2);fimg(end-10,end-10,end-10)=bar(2);end

fimg(10,10,10)=bar(1);
fimg(20,20,20)=bar(2);


config = struct('CameraPosition',[-4 0 0],...
    'CameraUpVector',[0,0,1],...
    'CameraTarget',[0 0 0],...
    'CameraViewAngle',15,...
    'BackgroundColor',[0 0 0]+1,...
    'Renderer','VolumeRendering',... %MaximumIntensityProjection%
    'Colormap',gray(256),...
    'Alphamap',(0:255)'/555,...
    'Lighting',0,...
    'IsosurfaceColor',[1 1 1],...
    'Isovalue',0.5);
config.Colormap = parula(256);
%gdmap = [zeros(128,1),(127:-1:0)'/127,zero(128,1)];
%drmap = [(0:127)'/127,ones(128,1),ones(128,1)];
%defaultMap = [gdmap;drmap;];
gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
% gdmap = [40/255+(255-40)/255*(0:127)'/127,125/255+(255-125)/255*(0:127)'/127,184/255+(255-184)/255*(0:127)'/127];
% domap = [ones(128,1),(63.5:0.5:127)'/127,(0:127)'/127];
defaultMap = [gdmap;flipud(drmap);];
% defaultMap = MY_viridis_colormap(256)/256;
config.Colormap = defaultMap;
config.Colormap = jet(256);
%config.Colormap = [ones(30,3)/1.3;jet(256-30)];
% config.Colormap = flipud([ones(256,1),(0:255)'/255,(0:255)'/255]);

% dorsal view
config.CameraPosition=[0 -3 0];
config.CameraUpVector = [-1 0 0];
F1=figure;volshow(fimg,config);
f1_rgb = getframe(F1);

% leteral view
config.CameraPosition=[0 +4 0];
config.CameraUpVector = [-1 0 0];
F2=figure;volshow(fimg,config);
f2_rgb = getframe(F2);

% anterior view
config.CameraPosition=[-3.45,-3.00,4.21]/2;
config.CameraUpVector = [-1 0 0];
F3=figure;volshow(fimg,config);
f3_rgb = getframe(F3);


config.CameraPosition=[-3 3 3]/1.6; %
config.CameraUpVector = [-1 0 0];%[0 0 1]
F4=figure;volshow(fimg,config);
f4_rgb = getframe(F4);


% hold on;volshow(timg,config);
% timg_rgb = getframe(F1);

config.CameraPosition=[-3 3 3]/1.6; %
config.CameraUpVector = [-1 0 0];%[0 0 1]
fimg(L<558 | L>1101)=nan;
fimgq = fimg;
fimgq(L<=690)=nan;
fimgq(10,10,10)=bar(1);
fimgq(20,20,20)=bar(2);
F3=figure;volshow(fimgq,config);
f3_rgb = getframe(F3);

end



