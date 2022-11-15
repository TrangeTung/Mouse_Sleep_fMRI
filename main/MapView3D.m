
clc;clear
codepath = 'F:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'F:\ECoG\dcc_yyl\';


path = 'F:\ECoG\Activation\';

Func3D = fullfile(path,'2nd','Xnrem_.nii');
Temp3D = fullfile(path,'2nd','Template_Mouse_v38.nii');
[f1_rgb,f2_rgb,f3_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,[0 5],codepath,1);

I = f1_rgb.cdata;%cat(1,f2_rgb.cdata,f3_rgb.cdata,f1_rgb.cdata);
for i=1:3
   I(:,:,i) = medfilt2(I(:,:,i),[5 5]);
end
imwrite(uint8(I),fullfile(path,'2nd','nrem.tiff'));



Func3D = fullfile(path,'2nd','Xrem_.nii');
Temp3D = fullfile(path,'2nd','Template_Mouse_v38.nii');
[f1_rgb,f2_rgb,f3_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,[0 30],codepath,1);

I = f1_rgb.cdata;%cat(1,f2_rgb.cdata,f3_rgb.cdata,f1_rgb.cdata);
for i=1:3
   I(:,:,i) = medfilt2(I(:,:,i),[5 5]);
end
imwrite(uint8(I),fullfile(path,'2nd','rem.tiff'));




Func3D = fullfile(path,'2nd','Xaw.nii');
Temp3D = fullfile(path,'Template_Mouse_v38.nii');
[f1_rgb,f2_rgb,f3_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,[-10 -1.75 1.75 10]/10,codepath,1);

I = f2_rgb.cdata;%cat(1,f2_rgb.cdata,f3_rgb.cdata,f1_rgb.cdata);
for i=1:3
   I(:,:,i) = medfilt2(I(:,:,i),[5 5]);
end
imwrite(uint8(I),fullfile(path,'2nd','AW.tiff'));











template = fullfile(path,'2nd','mean_T2.nii');
spmT_file = fullfile(path,'2nd','rem_raw.nii');
Colormap(   'statfile',spmT_file,...
            'bgmfile',template,...
            'slice',13:2:25,...
            'bar_value',[-40 -2.81 2.81 40],...
            'dest',fullfile(path,'2nd'),...
            'mapname','tmap_REM',...
            'cluster',10);
template = fullfile(path,'2nd','mean_T2.nii');
spmT_file = fullfile(path,'2nd','nrem_raw.nii');
Colormap(   'statfile',spmT_file,...
            'bgmfile',template,...
            'slice',13:2:25,...
            'bar_value',[-10 -1.81 4.81 20],...
            'dest',fullfile(path,'2nd'),...
            'mapname','tmap_NREM',...
            'cluster',10);
template = fullfile(path,'2nd','mean_T2.nii');
spmT_file = fullfile(path,'2nd','aw_raw.nii');
Colormap(   'statfile',spmT_file,...
            'bgmfile',template,...
            'slice',13:2:25,...
            'bar_value',[-20 -2.81 2.81 20],...
            'dest',fullfile(path,'2nd'),...
            'mapname','tmap_AW',...
            'cluster',10);













State = {'aw';'nrem';'rem'};
clear Rx
for sl=1:3
    k=0; img=[];
for il=1:46
    cd(fullfile(path,'1st',State{sl},num2str(il,'%02d')));
    spmT_dir = dir('spmT*.nii');
    if ~isempty(spmT_dir)
        k=k+1;
        img(:,:,:,k) = spm_read_vols(spm_vol(spmT_dir(1).name));
    end    
end
    Q = ~(mean(img,4)==0);
    V = fmask(img,Q);
    
    Cs = [5 10 20 30 size(V,2)];
    for cl=1:numel(Cs)
        [Rx(cl,sl)] = ICC(V(:,1:Cs(cl)), '1-k');
    end
end

F = figure;
plot(Rx,'-o');ylim([0 1])
xlim([0 6])

[X_,Y_,Z_] = sphere;
F = figure;
plot3([1,1,1,1,1]*+1,[1,2,3,4,5],Rx(:,1)*6); 
hold on;
for cl=1:numel(Cs);surf(X_*Rx(cl,1)/3+1,Y_*Rx(cl,1)/3+cl,Z_*Rx(cl,1)/3+Rx(cl,1)*6,'FaceColor','r','Edgecolor','none');end
plot3([1,1,1,1,1]*+3,[1,2,3,4,5],Rx(:,2)*6);
for cl=1:numel(Cs);surf(X_*Rx(cl,2)/3+3,Y_*Rx(cl,2)/3+cl,Z_*Rx(cl,2)/3+Rx(cl,2)*6,'FaceColor','g','Edgecolor','none');end
plot3([1,1,1,1,1]*+5,[1,2,3,4,5],Rx(:,3)*6);
for cl=1:numel(Cs);surf(X_*Rx(cl,3)/3+5,Y_*Rx(cl,3)/3+cl,Z_*Rx(cl,3)/3+Rx(cl,3)*6,'FaceColor','b','Edgecolor','none');end
grid on;
xlim([0 6]);zlim([0 6]);ylim([0 6])
axis equal
view([90 0]);


1;











