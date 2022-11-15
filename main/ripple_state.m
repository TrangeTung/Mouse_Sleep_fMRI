clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
cd(codepath)
WholePath = 'H:\ECoG\ripple-average\';

%% Quantified Event BOLD
close all
Excel0 = fullfile('H:\ECoG\Submit\v5','Supplementary_File1_Event_Quantitative_BOLD.xlsx');
Type = {'Awake ripple';'NREM ripple';...
    'spindle-coupled ripple';'spindle-uncoupled ripple';...
    'ripple-uncoupled spindle';'NREM spindle'};
MyMap = jet(256);
MyMap(128:129,:)=[0.7 0.7 0.7;0.7 0.7 0.7;];
F = figure('Position', [114 42 1731 953]);
for tl=1:numel(Type)
   [Data,~,~] = xlsread(Excel0, Type{tl});
    X = Data(2:end,:);
    X(isnan(X))=0;
    
    subplot(1,numel(Type),tl);
    imagesc(X);
    caxis([-0.3 0.3])
    colormap(MyMap);
    axis off
    title(Type{tl});
end
saveas(F,fullfile('H:\ECoG\Submit\v5\data','Supplementary_File1_Event_Quantitative_BOLD.emf'));

%% Time Series
close all
Excel0 = fullfile(WholePath,'timeseries.xlsx');
Regions = {'mPFC';'HPF';'Thal'};
Type = {'awake';'nrem';...
    'nrem_spindle_coupled_ripple';'nrem_spindle_uncoupled_ripple';...
    'nrem_ripple_uncoupled_spindle';'nrem_spindle'};
F = figure('Position', [169 51 1000 927]);
for rl=1:numel(Regions)
    [~,~,CellData] = xlsread(Excel0, Regions{rl});
    time=-12:2:12;
   for loop=1:3 
       subplot(3,3,(rl-1)*3+loop);
       switch loop
           case 1
               N1 = cell2mat(CellData(2:17,2:14));
               errorbar(-12:2:12,mean(N1,1),std(N1,0,1)/sqrt(size(N1,1)),...
                   'color',[239 126 027]/255,'linewidth',2,'CapSize',6);
               hold on;
               N2 = cell2mat(CellData(19:34,2:14));
               errorbar(-12:2:12,mean(N2,1),std(N2,0,1)/sqrt(size(N2,1)),...
                   'color',[020 129 055]/255,'linewidth',2,'CapSize',6);
               
               hv=zeros(size(N1,2),1);
               for tl=1:size(N1,2)
                  [hv(tl),~]=ttest(N1(:,tl),N2(:,tl)); 
               end
               plot(time(hv==1),10*ones(numel(find(hv==1)),1),'r*')
               hold off;
           case 2
               N1 = cell2mat(CellData(36:51,2:14));
               errorbar(-12:2:12,mean(N1,1),std(N1,0,1)/sqrt(size(N1,1)),...
                   'color',[153 000 204]/255,'linewidth',2,'CapSize',6);
               hold on;
               N2 = cell2mat(CellData(53:68,2:14));
               errorbar(-12:2:12,mean(N2,1),std(N2,0,1)/sqrt(size(N2,1)),...
                   'color',[204 153 255]/255,'linewidth',2,'CapSize',6);
               
               hv=zeros(size(N1,2),1);
               for tl=1:size(N1,2)
                  [hv(tl),~]=ttest(N1(:,tl),N2(:,tl)); 
               end
               plot(time(hv==1),10*ones(numel(find(hv==1)),1),'r*')
               hold off;
           case 3
               N1 = cell2mat(CellData(70:85,2:14));
               errorbar(-12:2:12,mean(N1,1),std(N1,0,1)/sqrt(size(N1,1)),...
                   'color',[173 127 069]/255,'linewidth',2,'CapSize',6);
               hold on;
               N1 = cell2mat(CellData(36:51,2:14));
               errorbar(-12:2:12,mean(N1,1),std(N1,0,1)/sqrt(size(N1,1)),...
                   'color',[153 000 204]/255,'linewidth',2,'CapSize',6);
               N1 = cell2mat(CellData(87:102,2:14));
               errorbar(-12:2:12,mean(N1,1),std(N1,0,1)/sqrt(size(N1,1)),...
                   'color',[243 024 016]/255,'linewidth',2,'CapSize',6);
               hold off;
       end
       ylim([-4 6]/10);xlim([-14.5 14.5])
       set(gca,'xtick',-12:6:12,'ticklength',[0.04 1],'linewidth',1.5);
       set(gca,'box','off','tickdir','out','fontsize',15)
   end
end
saveas(F,fullfile(WholePath,'TimeSeries.emf'));

%% Ripple pattern distance
V = spm_read_vols_4D(spm_vol(fullfile(codepath,'mean_epi.nii')));
lmask = V>800;
Func3D = fullfile(WholePath,'awake_ripple','awake_ripple.nii');
fMRI_4D = spm_read_vols_4D(spm_vol(Func3D));
V1 = fmask(fMRI_4D,lmask);V1(isnan(V1))=0; V1 = mean(V1,1);
Func3D = fullfile(WholePath,'nrem_ripple','nrem_ripple.nii');
fMRI_4D = spm_read_vols_4D(spm_vol(Func3D));
V2 = fmask(fMRI_4D,lmask);V2(isnan(V2))=0; V2 = mean(V2,1);
Func3D = fullfile(WholePath,'spindle_coupled','couple.nii');
fMRI_4D = spm_read_vols_4D(spm_vol(Func3D));
V3 = fmask(fMRI_4D,lmask);V3(isnan(V3))=0; V3 = mean(V3,1);
Func3D = fullfile(WholePath,'spindle_uncoupled','uncouple.nii');
fMRI_4D = spm_read_vols_4D(spm_vol(Func3D));
V4 = fmask(fMRI_4D,lmask);V4(isnan(V4))=0; V4 = mean(V4,1);
Func3D = fullfile(WholePath,'spindle','spindle.nii');
fMRI_4D = spm_read_vols_4D(spm_vol(Func3D));
V5 = fmask(fMRI_4D,lmask);V5(isnan(V5))=0; V5 = mean(V5,1);

[CC,px] = corr([V1(:),V2(:),V3(:),V4(:),V5(:)]);
figure;imagesc(CC);


%% Spatiotemporal motif
Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');
Func3D = fullfile(WholePath,'nrem_spindle','nrem_spindle.nii');
I = [];
for f=1:13
    Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                    'bgmfile',Temp3D,...
                    'slice',9:2:27,...
                    'bar_value',[-.1 -.02 .02 .1]*4,...
                    'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                    'dest',fullfile(WholePath,'spindle'),...
                    'mapname',['s'],...
                    'cluster',10);
    I = cat(1,I,Img_RGB);
end
imwrite(uint8(I),fullfile(WholePath,'nrem_spindle','nrem_spindle.tiff'));

Func3D = fullfile(WholePath,'spindle','spindle.nii');
I = [];
for f=1:13
    Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                    'bgmfile',Temp3D,...
                    'slice',9:2:27,...
                    'bar_value',[-.1 -.02 .02 .1]*4,...
                    'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                    'dest',fullfile(WholePath,'spindle'),...
                    'mapname',['s'],...
                    'cluster',10);
    I = cat(1,I,Img_RGB);
end
imwrite(uint8(I),fullfile(WholePath,'spindle','spindle.tiff'));


Func3D = fullfile(WholePath,'awake_ripple','awake_ripple.nii');
I = [];
for f=1:13
    Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                    'bgmfile',Temp3D,...
                    'slice',9:2:27,...
                    'bar_value',[-.1 -.02 .02 .1]*4,...
                    'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                    'dest',fullfile(WholePath,'awake_ripple'),...
                    'mapname',['s'],...
                    'cluster',10);
    I = cat(1,I,Img_RGB);
end
imwrite(uint8(I),fullfile(WholePath,'awake_ripple','awake_ripple.tiff'));

Func3D = fullfile(WholePath,'nrem_ripple','nrem_ripple.nii');
I = [];
for f=1:13
    Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                    'bgmfile',Temp3D,...
                    'slice',9:2:27,...
                    'bar_value',[-.1 -.02 .02 .1]*4,...
                    'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                    'dest',fullfile(WholePath,'awake_ripple'),...
                    'mapname',['s'],...
                    'cluster',10);
    I = cat(1,I,Img_RGB);
end
imwrite(uint8(I),fullfile(WholePath,'nrem_ripple','nrem_ripple.tiff'));

Func3D = fullfile(WholePath,'spindle_coupled','couple.nii');
I = [];
for f=1:13
    Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                    'bgmfile',Temp3D,...
                    'slice',9:2:27,...
                    'bar_value',[-.1 -.02 .02 .1]*4,...
                    'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                    'dest',fullfile(WholePath,'awake_ripple'),...
                    'mapname',['s'],...
                    'cluster',10);
    I = cat(1,I,Img_RGB);
end
imwrite(uint8(I),fullfile(WholePath,'spindle_coupled','couple.tiff'));

Func3D = fullfile(WholePath,'spindle_uncoupled','uncouple.nii');
I = [];
for f=1:13
    Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                    'bgmfile',Temp3D,...
                    'slice',9:2:27,...
                    'bar_value',[-.1 -.02 .02 .1]*4,...
                    'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                    'dest',fullfile(WholePath,'awake_ripple'),...
                    'mapname',['s'],...
                    'cluster',10);
    I = cat(1,I,Img_RGB);
end
imwrite(uint8(I),fullfile(WholePath,'spindle_uncoupled','uncouple.tiff'));

%% nrem Spindle
%
Func3D = fullfile(WholePath,'nrem_spindle','nrem_spindle.nii');
barlim = [-1 1]*0.50;
Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');
for f=1:13
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
    
    Ix = f2_rgb.cdata(81:350,71:490,:);
    Iz = f3_rgb.cdata(81:350,101:520,:);
    Iy = f1_rgb.cdata(:,101:460,:);
    Iwx = f4_rgb.cdata(:,41:540,:);
    Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);

    I1 = cat(2,Iy,permute(Ix,[2 1 3]));
    I2 = 255*ones(size(Iz,1),size(I1,2),3);
    I2(:,1:size(Iz,2),:)=Iz;
    I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
    
    I = cat(1,I1,I2);
    
    I1 = double(f4_rgb.cdata);
    I2 = double(f3_rgb.cdata);
    
    I = zeros(720,760,3);
    I1(end-100+1:end,end-100+1:end,:)=0;
    I(1:size(I1,1),1:size(I1,2),:) = I(1:size(I1,1),1:size(I1,2),:) + double(I1);
    I2(1:100,1:250,:)=0;
    I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:) = I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:)+double(I2);
    I(1:end-size(I2,1)+1,size(I1,2)+1:end,:) = 255;
    I(size(I1,1)+1:end,1:end-size(I2,2)+1,:) = 255;

    I1 = double(f1_rgb.cdata);
    I2 = double(f4_rgb.cdata);
    I3 = double(f3_rgb.cdata);
    I = cat(1,cat(1,I1,I2),I3);
    I = cat(1,I2,I3);
    dest = fullfile(WholePath,'nrem_spindle');
    mkdir(dest);
    
    imwrite(uint8(I),fullfile(dest,['Time_',num2str(f,'%02d'),'.tiff']));
    close all;
    

end

%% Spindle
%
Func3D = fullfile(WholePath,'spindle','spindle.nii');

Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');
for f=1:13
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
    
    Ix = f2_rgb.cdata(81:350,71:490,:);
    Iz = f3_rgb.cdata(81:350,101:520,:);
    Iy = f1_rgb.cdata(:,101:460,:);
    Iwx = f4_rgb.cdata(:,41:540,:);
    Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);

    I1 = cat(2,Iy,permute(Ix,[2 1 3]));
    I2 = 255*ones(size(Iz,1),size(I1,2),3);
    I2(:,1:size(Iz,2),:)=Iz;
    I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
    
    I = cat(1,I1,I2);
    
    I1 = double(f4_rgb.cdata);
    I2 = double(f3_rgb.cdata);
    
    I = zeros(720,760,3);
    I1(end-100+1:end,end-100+1:end,:)=0;
    I(1:size(I1,1),1:size(I1,2),:) = I(1:size(I1,1),1:size(I1,2),:) + double(I1);
    I2(1:100,1:250,:)=0;
    I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:) = I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:)+double(I2);
    I(1:end-size(I2,1)+1,size(I1,2)+1:end,:) = 255;
    I(size(I1,1)+1:end,1:end-size(I2,2)+1,:) = 255;

    I1 = double(f1_rgb.cdata);
    I2 = double(f4_rgb.cdata);
    I3 = double(f3_rgb.cdata);
    I = cat(1,cat(1,I1,I2),I3);
    I = cat(1,I2,I3);
    dest = fullfile(WholePath,'spindle');
    mkdir(dest);
    
    imwrite(uint8(I),fullfile(dest,['Time_',num2str(f,'%02d'),'.tiff']));
    close all;
    
%     Colormap(   'statfile',[Func3D,',',num2str(f)],...
%                 'bgmfile',Temp3D,...
%                 'slice',13:2:25,...
%                 'bar_value',[-.8 -.2 .2 .8],...
%                 'dest',dest,...
%                 'denoi_profile',fullfile(WholePath,'mask_RAM.nii'),...
%                 'mapname',['sPC_',num2str(f,'%02d')],...
%                 'cluster',10);
%
%     slice = 25:-2:13;
%
%     for sl=1:numel(slice)
%         I = spm_read_vols(spm_vol([Func3D,',',num2str(f)]));
%         J0 = MY_display_function_map_3D_nothreshold(I,[-.8 -.2 .2 .8],slice(sl));
%         Lmask = spm_read_vols(spm_vol(fullfile(WholePath,'mask_RAM.nii')));
%         L = flipud(Lmask(:,:,slice(sl))');
%         alpha=ones(size(L));alpha(L~=1)=0;
%         imwrite(uint8(J0),fullfile(dest,['slice_',num2str(slice(sl)),'.png']),'Alpha',alpha);
%
%     end
%

end
%}

%% awake ripple
%
Func3D = fullfile(WholePath,'awake_ripple','awake_ripple.nii');

Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');
for f=1:13
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
    
    Ix = f2_rgb.cdata(81:350,71:490,:);
    Iz = f3_rgb.cdata(81:350,101:520,:);
    Iy = f1_rgb.cdata(:,101:460,:);
    Iwx = f4_rgb.cdata(:,41:540,:);
    Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);

    I1 = cat(2,Iy,permute(Ix,[2 1 3]));
    I2 = 255*ones(size(Iz,1),size(I1,2),3);
    I2(:,1:size(Iz,2),:)=Iz;
    I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
    
    I = cat(1,I1,I2);
    
    I1 = double(f4_rgb.cdata);
    I2 = double(f3_rgb.cdata);
    
    I = zeros(720,760,3);
    I1(end-100+1:end,end-100+1:end,:)=0;
    I(1:size(I1,1),1:size(I1,2),:) = I(1:size(I1,1),1:size(I1,2),:) + double(I1);
    I2(1:100,1:250,:)=0;
    I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:) = I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:)+double(I2);
    I(1:end-size(I2,1)+1,size(I1,2)+1:end,:) = 255;
    I(size(I1,1)+1:end,1:end-size(I2,2)+1,:) = 255;

    I1 = double(f1_rgb.cdata);
    I2 = double(f4_rgb.cdata);
    I3 = double(f3_rgb.cdata);
    I = cat(1,cat(1,I1,I2),I3);
    I = cat(1,I2,I3);
    dest = fullfile(WholePath,'awake_ripple');
    mkdir(dest);
    
    imwrite(uint8(I),fullfile(dest,['Time_',num2str(f,'%02d'),'.tiff']));
    close all;
    
%     Colormap(   'statfile',[Func3D,',',num2str(f)],...
%                 'bgmfile',Temp3D,...
%                 'slice',13:2:25,...
%                 'bar_value',[-.8 -.2 .2 .8],...
%                 'dest',dest,...
%                 'denoi_profile',fullfile(WholePath,'mask_RAM.nii'),...
%                 'mapname',['sPC_',num2str(f,'%02d')],...
%                 'cluster',10);
%
%     slice = 25:-2:13;
%
%     for sl=1:numel(slice)
%         I = spm_read_vols(spm_vol([Func3D,',',num2str(f)]));
%         J0 = MY_display_function_map_3D_nothreshold(I,[-.8 -.2 .2 .8],slice(sl));
%         Lmask = spm_read_vols(spm_vol(fullfile(WholePath,'mask_RAM.nii')));
%         L = flipud(Lmask(:,:,slice(sl))');
%         alpha=ones(size(L));alpha(L~=1)=0;
%         imwrite(uint8(J0),fullfile(dest,['slice_',num2str(slice(sl)),'.png']),'Alpha',alpha);
%
%     end
%

end
%}

%% nrem ripple
%
Func3D = fullfile(WholePath,'nrem_ripple','nrem_ripple.nii');
Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');
for f=1:13
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
    
    Ix = f2_rgb.cdata(81:350,71:490,:);
    Iz = f3_rgb.cdata(81:350,101:520,:);
    Iy = f1_rgb.cdata(:,101:460,:);
    Iwx = f4_rgb.cdata(:,41:540,:);
    Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);

    I1 = cat(2,Iy,permute(Ix,[2 1 3]));
    I2 = 255*ones(size(Iz,1),size(I1,2),3);
    I2(:,1:size(Iz,2),:)=Iz;
    I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
    
    I = cat(1,I1,I2);
    
    I1 = double(f4_rgb.cdata);
    I2 = double(f3_rgb.cdata);
    
    I = zeros(720,760,3);
    I1(end-100+1:end,end-100+1:end,:)=0;
    I(1:size(I1,1),1:size(I1,2),:) = I(1:size(I1,1),1:size(I1,2),:) + double(I1);
    I2(1:100,1:250,:)=0;
    I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:) = I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:)+double(I2);
    I(1:end-size(I2,1)+1,size(I1,2)+1:end,:) = 255;
    I(size(I1,1)+1:end,1:end-size(I2,2)+1,:) = 255;

    I1 = double(f1_rgb.cdata);
    I2 = double(f4_rgb.cdata);
    I3 = double(f3_rgb.cdata);
    I = cat(1,cat(1,I1,I2),I3);
    I = cat(1,I2,I3);
    dest = fullfile(WholePath,'nrem_ripple');
    mkdir(dest);
    
    imwrite(uint8(I),fullfile(dest,['Time_',num2str(f,'%02d'),'.tiff']));
    close all;
    
%     Colormap(   'statfile',[Func3D,',',num2str(f)],...
%                 'bgmfile',Temp3D,...
%                 'slice',13:2:25,...
%                 'bar_value',[-.8 -.2 .2 .8],...
%                 'dest',dest,...
%                 'denoi_profile',fullfile(WholePath,'mask_RAM.nii'),...
%                 'mapname',['sPC_',num2str(f,'%02d')],...
%                 'cluster',10);
%
%     slice = 25:-2:13;
%
%     for sl=1:numel(slice)
%         I = spm_read_vols(spm_vol([Func3D,',',num2str(f)]));
%         J0 = MY_display_function_map_3D_nothreshold(I,[-.8 -.2 .2 .8],slice(sl));
%         Lmask = spm_read_vols(spm_vol(fullfile(WholePath,'mask_RAM.nii')));
%         L = flipud(Lmask(:,:,slice(sl))');
%         alpha=ones(size(L));alpha(L~=1)=0;
%         imwrite(uint8(J0),fullfile(dest,['slice_',num2str(slice(sl)),'.png']),'Alpha',alpha);
%
%     end
%

end


%}

%% spindle coupled ripple
%
Func3D = fullfile(WholePath,'spindle_coupled','couple.nii');

Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');
for f=1:13
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
    
    Ix = f2_rgb.cdata(81:350,71:490,:);
    Iz = f3_rgb.cdata(81:350,101:520,:);
    Iy = f1_rgb.cdata(:,101:460,:);
    Iwx = f4_rgb.cdata(:,41:540,:);
    Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);

    I1 = cat(2,Iy,permute(Ix,[2 1 3]));
    I2 = 255*ones(size(Iz,1),size(I1,2),3);
    I2(:,1:size(Iz,2),:)=Iz;
    I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
    
    I = cat(1,I1,I2);
    
    I1 = double(f4_rgb.cdata);
    I2 = double(f3_rgb.cdata);
    
    I = zeros(720,760,3);
    I1(end-100+1:end,end-100+1:end,:)=0;
    I(1:size(I1,1),1:size(I1,2),:) = I(1:size(I1,1),1:size(I1,2),:) + double(I1);
    I2(1:100,1:250,:)=0;
    I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:) = I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:)+double(I2);
    I(1:end-size(I2,1)+1,size(I1,2)+1:end,:) = 255;
    I(size(I1,1)+1:end,1:end-size(I2,2)+1,:) = 255;

    I1 = double(f1_rgb.cdata);
    I2 = double(f4_rgb.cdata);
    I3 = double(f3_rgb.cdata);
    I = cat(1,cat(1,I1,I2),I3);
    I = cat(1,I2,I3);
    dest = fullfile(WholePath,'spindle_coupled');
    mkdir(dest);
    
    imwrite(uint8(I),fullfile(dest,['Time_',num2str(f,'%02d'),'.tiff']));
    close all;
  

end
%}



%% spindle uncoupled ripple
%
Func3D = fullfile(WholePath,'spindle_uncoupled','uncouple.nii');

Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');
for f=1:13
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
    
    Ix = f2_rgb.cdata(81:350,71:490,:);
    Iz = f3_rgb.cdata(81:350,101:520,:);
    Iy = f1_rgb.cdata(:,101:460,:);
    Iwx = f4_rgb.cdata(:,41:540,:);
    Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);

    I1 = cat(2,Iy,permute(Ix,[2 1 3]));
    I2 = 255*ones(size(Iz,1),size(I1,2),3);
    I2(:,1:size(Iz,2),:)=Iz;
    I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
    
    I = cat(1,I1,I2);
    
    I1 = double(f4_rgb.cdata);
    I2 = double(f3_rgb.cdata);
    
    I = zeros(720,760,3);
    I1(end-100+1:end,end-100+1:end,:)=0;
    I(1:size(I1,1),1:size(I1,2),:) = I(1:size(I1,1),1:size(I1,2),:) + double(I1);
    I2(1:100,1:250,:)=0;
    I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:) = I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:)+double(I2);
    I(1:end-size(I2,1)+1,size(I1,2)+1:end,:) = 255;
    I(size(I1,1)+1:end,1:end-size(I2,2)+1,:) = 255;

    I1 = double(f1_rgb.cdata);
    I2 = double(f4_rgb.cdata);
    I3 = double(f3_rgb.cdata);
    I = cat(1,cat(1,I1,I2),I3);
    I = cat(1,I2,I3);
    dest = fullfile(WholePath,'spindle_uncoupled');
    mkdir(dest);
    
    imwrite(uint8(I),fullfile(dest,['Time_',num2str(f,'%02d'),'.tiff']));
    close all;
  

end
%}



