clc;clear
cd('H:\LFP_HIPP\hipp');
K = dir('m*');
g=28;%37,41,43,28
for ix=1:numel(K)
    cd(fullfile(K(ix).folder,K(ix).name));
    load(['s',num2str(g),'.mat']);
    load('fsEcg.mat');
    load('fsWav.mat');
    ch7=cell2mat(struct2cell(load('ch7.mat')));
    data = interp1((1:numel(ch7)) / 24414,ch7,1/1024:1/1024:numel(ch7) / 24414 , 'linear');
    %%   
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
    
%     x = (1:numel(y) )/1024;
% %     figure;plot(x,y_);hold on;plot(x,ye);
%     n=zeros(1,14745600);n(1,:)=threshold;
%     uC=zeros(1,14745600);uC(1,:)=mean(Cs);
%     uR=zeros(1,14745600);uR(1,:)=mean(Rs);
%     figure;plot(data2(1,9911050:9911357.2)/20+100);hold on;plot(y_(1,9911050:9911357.2));plot(ye(1,9911050:9911357.2));plot(n(9911050:9911357.2),'--');plot(ripple(9911050:9911357.2)*30+100);axis tight;
%     figure;plot(data2/20+100);hold on;plot(y_);plot(ye);plot(n,'--');plot(all_peak*30+100);axis tight;
%     
%     for i=1:14400
%         Cs1(((i-1)*1024+1):((i-1)*1024+1024),1)=Cs(i,1);
%         Rs1(((i-1)*1024+1):((i-1)*1024+1024),1)=Rs(i,1);
%     end
% 
%     figure;plot(data2(1,6288000:6298240)/10+200);hold on;plot(y_(1,6288000:6298240),'k');plot(ye(1,6288000:6298240),'w');plot(n(1,6288000:6298240),'--');plot(all_peak(1,6288000:6298240)*30+350);axis tight;%plot(Cs1(7898112:7908352,1)*100-200,'k');hold on;plot(Rs1(7898112:7908352,1)*20-150,'k');plot(uC(1,7898112:7908352)*100-200,'y--');plot(uR(1,7898112:7908352)*20-150,'y--');
%     figure;plot(Cs1(7898112:7908352,1)*2,'k');hold on;plot(Rs1(7898112:7908352,1)+1,'k');plot(uC(1,7898112:7908352)*2,'--');plot(uR(1,7898112:7908352)+1,'--');axis tight;%plot(all_state1(7898112:7908352,1));
%     figure;plot(data2(1,6296000:6297024)/10+200);hold on;plot(y_(1,6296000:6297024),'k');plot(ye(1,6296000:6297024),'y');plot(n(1,6296000:6297024),'--');plot(all_peak(1,6296000:6297024)*30+350);axis tight;

    
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

    %% find all peaks & 去掉低于阈值的peak，剩余peak标为2
      
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
%% brainstate

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
    
%     figure;plot(awake)
%%
    save(['ripple_peak',num2str(g),'.mat'],'ripple','awake_peak','nrem_peak','all_cut','awake_cut','nrem_cut');
    
    g=g+1;
    
end