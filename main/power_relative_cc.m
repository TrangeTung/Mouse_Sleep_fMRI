%% spectrum
de_data=load('de_data.mat');
movingwin=[3 1]; % set the moving window dimensions
params.Fs=fs; % sampling frequency
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
%% relative change
for idx=[1:27]
    tt1=cell2mat(struct2cell(load(['I:\fig2_cc_validation\roi_3\roi_cortex_x_1009\rsp\roi',num2str(idx),'.mat'],'cortex_mean')));
    
    len=epoch(:,2)-epoch(:,1);
    ll=find(len<=15);
    epoch(ll,:)=[];
    
    load(['C:\Users\yuyalin\Desktop\result_jin\state',num2str(idx),'.mat']);
    SE = strel('disk',5);
    Bs = imerode(state_tr,SE);
    bold_base=mean(tt1(1,find(Bs==1)));
    
    for ii=1:size(epoch,1)
        test=tt1(epoch(ii,1):epoch(ii,2));
        nrem_bold(ii)=mean(test);
        clear test
    end
    nrem_bold1=(nrem_bold-bold_base)/bold_base;
    
    s=0;
    for k=[2:8 10:16]%2:8 10:16
        s1=cell2mat(struct2cell(load(['I:\ecog\x_S_100_',num2str(k),'.mat'],'S')));
        %figure;plot_matrix(s1(:,5:117),t,f(5:117));caxis([-10 20]);
        s=s+s1;
    end
    s=s./14;  %

    s1=s(:,1:21);  %1-4
    s2=s(:,22:33);  %5-10
    s3=s(:,34:57);   %11-20
    s4=s(:,57:157);  %21-40
    s5=s(:,158:337);  %41-100
    %load(['C:\Users\yuyalin\Desktop\result_jin\states_x\state_check\s',num2str(idx),'.mat']);%
    all_state=round(resample(state_tr,2,1));    SE = strel('disk',10); all_state = imerode(all_state,SE);
    ecog_base=s1(find(all_state==1),:);ecog_base1=mean(mean(ecog_base)); clear ecog_base
    ecog_base=s2(find(all_state==1),:);ecog_base2=mean(mean(ecog_base));clear ecog_base
    ecog_base=s3(find(all_state==1),:);ecog_base3=mean(mean(ecog_base));clear ecog_base
    ecog_base=s4(find(all_state==1),:);ecog_base4=mean(mean(ecog_base));clear ecog_base
    ecog_base=s5(find(all_state==1),:);ecog_base5=mean(mean(ecog_base));clear ecog_base
    for ii=1:size(epoch1,1)
        test1=s1(epoch1(ii,1):epoch1(ii,2),:);nrem_delta(ii)=mean(mean(test1)); 
        test2=s2(epoch1(ii,1):epoch1(ii,2),:);nrem_theta(ii)=mean(mean(test2)); 
        test3=s3(epoch1(ii,1):epoch1(ii,2),:);nrem_alpha(ii)=mean(mean(test3)); 
        test4=s4(epoch1(ii,1):epoch1(ii,2),:);nrem_beta(ii)=mean(mean(test4)); 
        test5=s5(epoch1(ii,1):epoch1(ii,2),:);nrem_gamma(ii)=mean(mean(test5)); 
    end
    nrem_delta1=nrem_delta-ecog_base1;
    nrem_theta1=nrem_theta-ecog_base2;
    nrem_alpha1=nrem_alpha-ecog_base3;
    nrem_beta1=nrem_beta-ecog_base4;
    nrem_gamma1=nrem_gamma-ecog_base5;
    save(['I:\ecog_re-ref\rela\nrem\delta\cortex',num2str(idx),'.mat'],'nrem_delta1','nrem_bold1');
    save(['I:\ecog_re-ref\rela\nrem\theta\cortex',num2str(idx),'.mat'],'nrem_theta1','nrem_bold1');
    save(['I:\ecog_re-ref\rela\nrem\alpha\cortex',num2str(idx),'.mat'],'nrem_alpha1','nrem_bold1');
    save(['I:\ecog_re-ref\rela\nrem\beta\cortex',num2str(idx),'.mat'],'nrem_beta1','nrem_bold1');
    save(['I:\ecog_re-ref\rela\nrem\gamma\cortex',num2str(idx),'.mat'],'nrem_gamma1','nrem_bold1');
    clearvars -except 'path' 'idx';
end
%%
bold=[];
ecog_gamma=[];
for ii=[1:9 13:15 17:19 22:27]
    load(['I:\ecog_re-ref\rela\rem\gamma\cortex',num2str(ii),'.mat']);
    bold=[bold,nrem_bold1];
    ecog_gamma=[ecog_gamma,nrem_gamma1];
    clear nrem_bold1
    clear nrem_ecog1
end
 %%
delta=[ecog_delta;bold];%a=find(delta(1,:)<0);delta(:,a)=[];clear a
theta=[ecog_theta;bold];%a=find(theta(1,:)<0);theta(:,a)=[];clear a
alpha=[ecog_alpha;bold];%a=find(alpha(1,:)<0);alpha(:,a)=[];clear a
beta=[ecog_beta;bold];%a=find(beta(1,:)<0);beta(:,a)=[];clear a
gamma=[ecog_gamma;bold];%a=find(gamma(1,:)<0);gamma(:,a)=[];clear a
%% cc
[R1,P1]=corrcoef(delta(1,:),delta(2,:));p1=polyfit(delta(1,:),delta(2,:),1);yfit1=polyval(p1,delta(1,:));
[R2,P2]=corrcoef(theta(1,:),theta(2,:));p2=polyfit(theta(1,:),theta(2,:),1);yfit2=polyval(p2,theta(1,:));
[R3,P3]=corrcoef(alpha(1,:),alpha(2,:));p3=polyfit(alpha(1,:),alpha(2,:),1);yfit3=polyval(p3,alpha(1,:));
[R4,P4]=corrcoef(beta(1,:),beta(2,:));p4=polyfit(beta(1,:),beta(2,:),1);yfit4=polyval(p4,beta(1,:));
[R5,P5]=corrcoef(gamma(1,:),gamma(2,:));p5=polyfit(gamma(1,:),gamma(2,:),1);yfit5=polyval(p5,gamma(1,:));
%%
figure;subplot(1,5,1);scatter(delta(1,:),delta(2,:));hold on;plot(delta(1,:),yfit1);xlabel('delta');ylabel('hipp bold');
subplot(1,5,2);scatter(theta(1,:),theta(2,:));hold on;plot(theta(1,:),yfit2);xlabel('theta');ylabel('hipp bold');
subplot(1,5,3);scatter(alpha(1,:),alpha(2,:));hold on;plot(alpha(1,:),yfit3);xlabel('alpha');ylabel('hipp bold');
subplot(1,5,4);scatter(beta(1,:),beta(2,:));hold on;plot(beta(1,:),yfit4);xlabel('beta');ylabel('hipp bold');
subplot(1,5,5);scatter(gamma(1,:),gamma(2,:));hold on;plot(gamma(1,:),yfit5);xlabel('gamma');ylabel('hipp bold');