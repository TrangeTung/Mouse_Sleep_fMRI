%% events_avg
for idx=[28:46]
fpath=[path{idx},'Results\',num2str(Animal_EPI_folder{idx}),'\snrrsm2dseq.nii'];
fmri=spm_vol(fpath);
fdata=spm_read_vols(fmri); %
mdata=spm_read_vols(spm_vol('D:\lianglab\code\fmri\mouse\lmask_Mouse_v38.nii'));
data=fmask(fdata,mdata);
clear fdata

%data1=(data-mean(data')')./std(data')';

% load(['C:\Users\yuyalin\Desktop\result_jin\states_x\state_tr\s',num2str(idx),'.mat']);
% load(['I:\events_0821\spindle\spindle_peak',num2str(idx),'.mat'],'cut');
% load(['C:\Users\yuyalin\Desktop\ripple_peak\ripple_peak',num2str(idx),'.mat'],'awake_cut','nrem_cut');
% 
% spin=cut./2048;spin1=unique(round(spin(:,2)));
% awake_ripp=awake_cut./2048;awake_ripp1=unique(round(awake_ripp(:,2)));
% nrem_ripp=nrem_cut./2048;nrem_ripp1=unique(round(nrem_ripp(:,2)));
% clear 'awake_cut' 'awake_ripp' 'cut' 'nrem_cut' 'nrem_ripp' 'spin'
%% nrem ripple
path1=['I:\events_0821\baseline1_x\nrem\',num2str(idx),'\'];
mkdir(path1);
load(['C:\Users\yuyalin\Desktop\ripple_peak\ripple_peak',num2str(idx),'.mat'],'nrem_cut');
n1=nrem_cut./2048;
n_ri=unique(round(n1(:,2)));
n_ri(find(n_ri<7))=[];
n_ri(find(n_ri>7194))=[];

base2=[];
for i=1:length(n_ri)
    base1=data(:,(n_ri(i)-4):(n_ri(i)-2));
    base2=[base2,base1];
    clear base1
end

base=mean(base2,2);
clear base2
navg=[];
for i=1:length(n_ri)
    navg(:,:,i)=(data(:,(n_ri(i)-6):(n_ri(i)+6))-base)./base;   
end
navg1=sum(navg,3);
ss=i;
save([path1,'all.mat'],'ss');

cc=cell(1,13);
for k=1:13
    cc{k}=funmask(navg1(:,k),mdata);
   % cc{k}=reshape(navg1(:,k),114,85,35);
end

fmri1=fmri(1);
for j=1:13
    fmri1.fname=[path1,'avg',num2str(j),'.nii']; %
    fmri1.dt=[16,0];
    spm_write_vol(fmri1,cc{j});
end
end
%% frequency_cc
k=1;
for ii=[28:46]
   for j=1:13
nrem_bold(:,:,:,j)=spm_read_vols(spm_vol(['I:\events_0821\state_avg1_qu_events\diff\nrem\',num2str(ii),'\avg',num2str(j),'.nii']));
awake_bold(:,:,:,j)=spm_read_vols(spm_vol(['I:\events_0821\state_avg1_qu_events\diff\awake\',num2str(ii),'\avg',num2str(j),'.nii']));
    end
    
% hh=load(['I:\events_pattern_0630\nrem_ripple\',num2str(ii),'\all.mat']);hh=hh.ss;
% hh1=load(['I:\events_pattern_0630\awake_ripple\',num2str(ii),'\all.mat']);hh1=hh1.ss;

hipp_mask=spm_read_vols(spm_vol('I:\events_validation\roi_mask\hipp_1.nii'));
mpfc_mask=spm_read_vols(spm_vol('I:\events_validation\roi_mask\mpfc_1.nii'));

mpfc_awake=mean(fmask(awake_bold,mpfc_mask));
mpfc_nrem=mean(fmask(nrem_bold,mpfc_mask));
mpfc_diff=mpfc_nrem-mpfc_awake;

hipp_awake=mean(fmask(awake_bold,hipp_mask));
hipp_nrem=mean(fmask(nrem_bold,hipp_mask));
hipp_diff=hipp_nrem-hipp_awake;
clearvars -except 'mpfc_diff' 'hipp_diff' 'ii' 'hrf' 'k' 'mpfc_cc' 'mpfc_p' 'hipp_cc' 'hipp_p';
%%
movingwin=[1 0.1]; % set the moving window dimensions
params.Fs=1024; % sampling frequency
params.fpass=[1 100]; % frequencies of interest
params.tapers=[3 5]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
load(['I:\events_validation\cross_cc_0815\power\mpfc\awake\mpfc',num2str(ii-27),'.mat']);  %%
%ss=ss(:,10*1024+1:18*1024);
for j=1:size(ss,1)
    [S3(:,:,j),t2,f2]=mtspecgramc(ss(j,:),movingwin,params);
end
awake_power=mean(S3,3);
clear ss S3;
load(['I:\events_validation\cross_cc_0815\power\mpfc\nrem\mpfc',num2str(ii-27),'.mat']);   %%
%ss=ss(:,10*1024+1:18*1024);
for j=1:size(ss,1)
    [S3(:,:,j),t2,f2]=mtspecgramc(ss(j,:),movingwin,params);
end
nrem_power=mean(S3,3);
clear ss S3;
 %%
awake1=10*log10(awake_power);
nrem1=10*log10(nrem_power);
diff1=nrem1-awake1;
for i=1:100
diff11(:,i)=conv(diff1(:,i),hrf);
end
diff11=diff11(1:251,:);

%diff1=diff1(:,4:end);
diff2=resample(diff11,13,251);%13,251 4,71
diff2=resample(diff2,260,13);

%hipp_diff=resample(hipp_diff,260,13);
mpfc_diff=resample(mpfc_diff,260,13);
% mpfc_cc=corr(diff2,mpfc_diff(6:9)');
% mpfc_cc1=corr(mean(diff2,2),mpfc_diff(6:9)');
% hipp_cc=corr(diff2,hipp_diff(6:9)');
% hipp_cc1=corr(mean(diff2,2),hipp_diff(6:9)');
%hipp_diff=resample(hipp_diff,260,13);


[mpfc_cc(:,k),mpfc_p(:,k)]=corr(diff2,mpfc_diff');
%[hipp_cc(:,k),hipp_p(:,k)]=corr(diff2,hipp_diff');
%mpfc_cc1(:,k)=corr(diff2,circshift(mpfc_diff',-20));
k=k+1;

% kk=1;
% for k=100:-1:-100
% hipp_cc(:,kk)=corr(diff2,circshift(hipp_diff',k));
% hipp_cc1(kk)=corr(mean(diff2,2),circshift(hipp_diff',k));
% kk=kk+1;
% end

%save(['I:\events_0821\diff\diff\mpfc',num2str(ii-27),'.mat'],'mpfc_cc','mpfc_cc1');
clearvars -except 'ii' 'hrf' 'k' 'mpfc_cc' 'mpfc_p' 'hipp_cc' 'hipp_p';
end
