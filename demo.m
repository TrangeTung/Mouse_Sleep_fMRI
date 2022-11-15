%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fistly
% Download code packages from the links ("download_URL.txt" in each folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% codepath = '....\..\rsfMRI\';
% filepath = '....\..\examplefile\';
addpath(genpath(codepath));
%% eletrophysiological denoising
tic;
filename = fullfile(filepath,'raw_data.xlsx');
ecg = xlsread(filename);
EEG.srate=fix([fullfile(filepath,'fs.xlsx')]);
EEG.pnts=length(ecg);
nSlice=22;
i=1;
for a=1:nSlice:size(startTrigger,1)
    trig(i)=startTrigger(a,1)/fsWav*fsEcg;
    i=i+1;
end
filename = fullfile(filepath,'raw_data.xlsx');
trig = xlsread(filename);
trig=ceil(trig);
clear data;
for i=[1:16]
    EEG.data=ecg(i,:);
    de_data=fmrib_fastr(EEG,0,1,30,trig,0,0,0);
    name_res=['ch',num2str(i+2),'.mat'];
    save(name_res,'de_data','-v7.3');
    clear de_data;
end
de_data=resample(de_data,1024,24414);
%% power spectrum
filename = fullfile(filepath,'denoised_ECoG_data.xlsx');
de_data = xlsread(filename);
movingwin=[3 1]; % set the moving window dimensions
params.Fs=1024; % sampling frequency
params.fpass=[1 30]; % frequencies of interest
params.tapers=[3 5]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
%data=ecg(dataStart:dataEnd,channel); % data from channel 
[S3,t2,f2]=mtspecgramc(de_data,movingwin,params);
figure;
plot_matrix(S3,t2,f2);
xlabel('Time (sec)');ylabel('Frequency(Hz)'); % plot spectrogram
colormap('jet');
caxis([0 40]); colorbar; 
%% event detection
filename = fullfile(filepath,'denoised_LFP_data.xlsx');
data = xlsread(filename);
movingwin=[1 1]; % set the moving window dimensions
params.Fs=1024; % sampling frequency
params.fpass=[120 350]; % frequencies of interest
params.tapers=[3 5]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
%data=ecg(dataStart:dataEnd,channel); % data from channel
[S3,t2,f2]=mtspecgramc(data,movingwin,params);
%     figure;
%     plot_matrix(S3,t2,f2);
%     xlabel('Time (sec)');ylabel('Frequency(Hz)'); % plot spectrogram
%     colormap('jet');
%     caxis([-10 20]); colorbar;

% for motion
Cs = mean(S3(:,f2>250),2);
Bs  = Cs>mean(Cs);%mean(Cs)
SE = strel('disk',1);
Bs = imdilate(Bs,SE);
% for ripple
Rs = mean(S3(:,f2>120&f2<250),2);
Rs(Bs) = 0.2*mean(Rs);%1
Rp = Rs> mean(Rs);%1
Rp2 = imdilate(Rp,SE);
Rp2E = interp1((1:numel(Rp2)),double(Rp2), 1/1024:1/1024:numel(Rp2) ,'pchip');
% for amplitude
[b,a] = butter(10,120/512,'high');
y_ = filtfilt(b,a,data2);
[b,a] = butter(10,250/512,'low');
y_ = filtfilt(b,a,y_);
ye = envelope(y_,100,'rms');
m=mean(ye);
threshold=mean(ye)+std(ye)*3;%3
RpA = ye>threshold;
NUM = min([numel(RpA),numel(Rp2E)]);
RpB = RpA(1:NUM).*Rp2E(1:NUM);
figure;
plot_matrix(S3,t2,f2);
xlabel('Time (sec)');ylabel('Frequency(Hz)'); % plot spectrogram
colormap('jet');
caxis([-10 20]); colorbar;
%     hold on;plot(t2,Bs*10+200,'w','linewidth',2)
%     plot(t2,Rp*10+170,'k','linewidth',2)
%     plot(x,y_/10+100,'k','linewidth',0.5);
hold on;plot((1:numel(all_peak) )/1024,all_peak*10+80,'k','linewidth',2);
hold on;plot((1:numel(RpB) )/1024,RpB*10+80,'k','linewidth',2)
ylim([50 250])

%     a=length(find(gradient(RpB)>0.1));

%% find all peaks & excluded beloe the threshold
[P,L]=findpeaks(ye);
peak=ones(1,14745600);
peak(L)=2;
for k=1:length(L)
    y=L(k);
    if ye(y) < threshold
        peak(y)=1;
    else
        k=k+1;
    end
end

all_peak=RpB.*peak;
all_peak(1,14399*1024:14400*1024)=0;
%     figure;plot(all_peak*100);hold on;plot(env);plot(rip1);
%% 以peak为中心，寻找起始点和终点，并删除
p=find(all_peak==2)';
if p(1,1)<m
    p(1,:)=[];
end
if p(1,1)<m
    p(1,:)=[];
end

%     if 14745600-p(end,1)<m
%         p(end,:)=[];
%     end
%     if 14745600-p(end,1)<m
%         p(end,:)=[];
%     end

ripple=zeros(1,14745600);
n=1;
for i=1:length(p)
    x=1;
    for x=1:p(i,1)
        bg=p(i,1)-x;
        if ye(1,bg)<m
            break
        else
            x=x+1;
        end
    end
    
    x=1;
    for x=1:p(i,1)
        ed=p(i,1)+x;
        if ye(1,ed)<m
            break
        else
            x=x+1;
        end
    end
    if ed-bg<3
        break
    else
        piece(n,1)=bg;
        piece(n,3)=ed;
        ripple(1,bg:ed)=1;
        n=n+1;
    end
end

%     figure;plot(ripple*30+30);hold on;plot(rip1);plot(env);

all_cut=unique(piece,'rows');

for i=1:length(all_cut)
    f1=all_cut(i,1);
    f2=all_cut(i,3);
    [Y,I]=max(y_(1,f1:f2));
    all_cut(i,2)=I+f1;
end

ripple(all_cut(:,2))=2;
%     figure;plot(ripple*30+30);hold on;plot(rip1/5);plot(env/5);
%% brainstate label

for i=1:14400
    all_state1(((i-1)*1024+1):((i-1)*1024+1024),1)=all_state(i,1);
end

e=zeros(14745600,1);
e1=zeros(14745600,1);
e(find(all_state1==1))=1;
e1(find(all_state1==3))=1;
awake_peak=ripple.*e';
nrem_peak=ripple.*e1';

awake_cut=[];
nrem_cut=[];
for i=1:length(all_cut)
    if awake_peak(all_cut(i,1))==1
        awake_cut=[awake_cut;all_cut(i,:)];
    else
        nrem_cut=[nrem_cut;all_cut(i,:)];
    end
end
%%
save(['ripple_peak',num2str(g),'.mat'],'ripple','awake_peak','nrem_peak','all_cut','awake_cut','nrem_cut');
g=g+1;

