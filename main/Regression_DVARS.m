clc;clear
WholePath = 'F:\ECoG\shang_m2_3';
filename = fullfile(WholePath,'snrsm2dseq.nii');
hdr = spm_vol(filename);
fMRI_4D = spm_read_vols_4D(hdr);

lmask = spm_read_vols_4D(spm_vol(fullfile(WholePath,'mask.nii')));

rp = load(fullfile(WholePath,'rp_sm2dseq.txt'));
drp = rp - [rp(1,:);rp(1:end-1,:)];
PCs = load(fullfile(WholePath,'PC40.txt'));


RegBasFunc1 = [ones([length(rp),1]) rp drp ];
RegBasFunc2 = [ones([length(rp),1]) rp drp PCs(:,1:10)];
RegBasFunc3 = [ones([length(rp),1]) rp drp PCs(:,1:40)];


data_ready_regress0 = fmask(fMRI_4D,lmask);   
data_ready_regress1 = data_ready_regress0*0;
data_ready_regress2 = data_ready_regress0*0;
data_ready_regress3 = data_ready_regress0*0;

for iii = 1:size(data_ready_regress0,1)
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc1);
    data_ready_regress1(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc2);
    data_ready_regress2(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc3);
    data_ready_regress3(iii,:) = Residual + Beta(1);
end

data = data_ready_regress0;
dd = data-[data(:,1),data(:,1:end-1)];
DVARS(:,1) =rms(dd);
data = data_ready_regress1;
dd = data-[data(:,1),data(:,1:end-1)];
DVARS(:,2) =rms(dd);
data = data_ready_regress2;
dd = data-[data(:,1),data(:,1:end-1)];
DVARS(:,3) = rms(dd);
data = data_ready_regress3;
dd = data-[data(:,1),data(:,1:end-1)];
DVARS(:,4) = rms(dd);

drp(1:3)=drp(1:3)/20;
drp(4:6)=drp(4:6)*5;
FD=sqrt(drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2);

X = corr(FD(:),data_ready_regress1');
[~,s] = sort(X,'descend');

F = figure('Position', [680 401 519 577]);
subplot(10,1,1);plot(FD*100);xlim([0 7200]);
subplot(10,1,2);plot(DVARS);xlim([0 7200]);
data=(data_ready_regress0-mean(data_ready_regress0,2))./std(data_ready_regress0,0,2);
subplot(10,1,3:4);
plot(nanmean(data,1));xlim([0 7200]);ylim([-2 2])
% imagesc(data(s,:));caxis([-2 2]);colormap('gray');
data=(data_ready_regress1-mean(data_ready_regress1,2))./std(data_ready_regress1,0,2);
subplot(10,1,5:6);
plot(nanmean(data,1));xlim([0 7200]);ylim([-2 2])
% imagesc(data(s,:));caxis([-2 2]);colormap('gray');
data=(data_ready_regress2-mean(data_ready_regress2,2))./std(data_ready_regress2,0,2);
subplot(10,1,7:8);
plot(nanmean(data,1));xlim([0 7200]);ylim([-2 2])
% imagesc(data(s,:));caxis([-2 2]);colormap('gray');
data=(data_ready_regress3-mean(data_ready_regress3,2))./std(data_ready_regress3,0,2);
subplot(10,1,9:10);
plot(nanmean(data,1));xlim([0 7200]);ylim([-2 2])
% imagesc(data(s,:));caxis([-2 2]);colormap('gray');

