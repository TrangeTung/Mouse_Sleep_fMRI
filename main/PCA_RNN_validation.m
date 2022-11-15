clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\raw_data_yyl';

ihdr = spm_vol(fullfile('H:\ECoG\pca\','mask_RAM.nii'));
lmask = spm_read_vols(ihdr);
ROI4Type = {'SSs';'AUD';'AMY';'VIS';'PnV'};
ROI5Type = {'vThal';'VMH';'VTA';'LC';};
ROI6Type = {'RSP';'ACA';'dThal';'TT';'PIR';'MnR';'ANcr2';'PF'};
ROInames = [ROI4Type;ROI5Type;ROI6Type];

%{
ALLmask = zeros([size(lmask),numel(ROInames)]);al=0;
ROIName = fullfile('D:\ECoG\timeseries_roi_213\PCA_TEST\Sensitivity','ROI_4.nii');
I = spm_read_vols(spm_vol(ROIName));
for it=1:numel(ROI4Type);al=al+1;Ix=I==it;Ix=Ix+flip(Ix,1);Ix(Ix~=0)=1;ALLmask(:,:,:,al)=Ix;end
ROIName = fullfile('D:\ECoG\timeseries_roi_213\PCA_TEST\Sensitivity','ROI_5.nii');
I = spm_read_vols(spm_vol(ROIName));
for it=1:numel(ROI5Type);al=al+1;Ix=I==it;Ix=Ix+flip(Ix,1);Ix(Ix~=0)=1;ALLmask(:,:,:,al)=Ix;end
ROIName = fullfile('D:\ECoG\timeseries_roi_213\PCA_TEST\Sensitivity','ROI_6.nii');
I = spm_read_vols(spm_vol(ROIName));
for it=1:numel(ROI6Type);al=al+1;Ix=I==it;Ix=Ix+flip(Ix,1);Ix(Ix~=0)=1;ALLmask(:,:,:,al)=Ix;end
Rv = fmask(ALLmask,lmask);
for idx = [11:30 1:10 31:46]
   dest = fullfile(WholePath, num2str(idx));
   filename = fullfile(dest,'sLnrrsm2dseq.nii');
   fMRI_4D = spm_read_vols_4D(spm_vol(filename));
   Sv = fmask(fMRI_4D,lmask); clear fMRI_4D
   S = (Sv'*Rv)';
   Sz = (S-mean(S,2))./std(S,0,2);
   save(fullfile(dest,'ROIseries.mat'),'S','Sz')
   clear Sv
end
%}

%
Signal = [];
State = [];
for idx=1:46
  
   dest = fullfile(WholePath, num2str(idx));
   load(fullfile(dest,'ROIseries.mat'))
   E = Sz'; 
   
    [b,a]=butter(3,0.1);
    for ix=1:size(E,2); F(:,ix) = filtfilt(b,a,E(:,ix)'); end
    
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));
    b1 = interp1(1:14400,all_state,2:2:14400,'nearest')';
    b1(b1==3)=2;
    b1(b1==4)=3;
    
    Signal = [Signal;E];
    State = [State;b1];

end

Cs = 60; % cut TR
s12_=0;s23_=0;s31_=0;s21_=0;
s11_=0;s22_=0;s33_=0;
for sl=Cs+1:numel(State)-Cs
    if State(sl)==1 & State(sl+1)==2
        if numel(unique(State(sl+(-Cs:-1))))==1
            s12_=s12_+1;XX12_(s12_)=sl;
        end
    end
    if State(sl)==2 & State(sl+1)==3
            s23_=s23_+1;XX23_(s23_)=sl;
    end
    if State(sl)==3 & State(sl+1)==1
            s31_=s31_+1;XX31_(s31_)=sl;
    end
    if State(sl)==2 & State(sl+1)==1
        if numel(unique(State(sl+(-Cs:-1))))==1
            s21_=s21_+1;XX21_(s21_)=sl;
        end
    end
    
    y = State(sl+(-Cs:3));
    y_ = unique(y);
    if numel(y_)==1;
        if y_==1; s11_=s11_+1;XX11_(s11_)=sl;end
        if y_==2; s22_=s22_+1;XX22_(s22_)=sl;end
        if y_==3; s33_=s33_+1;XX33_(s33_)=sl;end
    end
    
end

clear Y*_ S 
delete(fullfile(WholePath,'WO*.mat'));
Type = {'11';'22';'33';'12';'21';'23';'31'};
for tl=1:numel(Type)
    bins = 5;
    if tl==3 | tl==7| tl==6;bins=25;end
    Y=[];S=[];
    eval(['sx=s',Type{tl},'_;']);
    eval(['XX=XX',Type{tl},'_;']);
    
    BINS = [13642;9000;12560;713;236;167;165];
    
    if tl==1;StateIndx = 1*ones(Cs*2+1,1);end
    if tl==2;StateIndx = 2*ones(Cs*2+1,1);end
    if tl==3;StateIndx = 3*ones(Cs*2+1,1);end
    
    if tl==4;StateIndx = [1*ones(Cs+1,1);2*ones(Cs,1)];end
    if tl==5;StateIndx = [2*ones(Cs+1,1);1*ones(Cs,1)];end
    if tl==6;StateIndx = [2*ones(Cs+1,1);3*ones(Cs,1)];end
    if tl==7;StateIndx = [3*ones(Cs+1,1);1*ones(Cs,1)];end
    
    if tl==4;StateIndx = [1*ones(Cs+1,1);4*ones(Cs,1)];end
    if tl==5;StateIndx = [2*ones(Cs+1,1);5*ones(Cs,1)];end
    if tl==6;StateIndx = [2*ones(Cs+1,1);6*ones(Cs,1)];end
    if tl==7;StateIndx = [3*ones(Cs+1,1);7*ones(Cs,1)];end
    
    for loop=1:BINS(tl)
        M=[]; 
        XX = XX(randperm(numel(XX)));
        for sl=1:bins
            M(:,:,sl) = Signal(XX(sl)+(-Cs:Cs),:);
        end
       Y(:,:,loop) = mean(M,3);
       
       S(1,:,loop) = StateIndx;
    end
    
    eval(['Y',Type{tl},'_=Y;']);
    
    dest = fullfile(WholePath);
    if ~exist(dest,'dir');mkdir(dest);end
    
    save(fullfile(dest,['WO_Signal',Type{tl},'.mat']),'Y');
    save(fullfile(dest,['WO_State',Type{tl},'.mat']),'S');
end

%}


Type = {'11';'22';'33';'12';'21';'23';'31'};

ColorX={[228,147,18]/255;[17,140,56]/255;[9,102,164]/255};
for bs=4:7
    
    if bs==4;Pks=[4 1];ColorVec=ColorX([2,1]);end
    if bs==5;Pks=[5 2 6];ColorVec=ColorX([1,2,3]);end % 6 3
    if bs==6;Pks=[6 5 2];ColorVec=ColorX([3,1,2]);end
    if bs==7;Pks=[7 3];ColorVec=ColorX([1,3]);end
    
    clear x* y* Signal State


    clear x* y* Signal State
    for tl = 1:numel(Pks)
        load(fullfile(WholePath,['WO_Signal',Type{Pks(tl)},'.mat']));
        Signal = Y;
        load(fullfile(WholePath,['WO_State',Type{Pks(tl)},'.mat']));
        State = S;
        N = Pks(tl)*ones(size(S,3),1);
        unique(N)
        eval(['x',num2str(tl),'=Signal;']);
    end

    close all
    
    for ps=1:size(x1,2)
        
        y1=squeeze(x1(:,ps,:));
        y2=squeeze(x2(:,ps,:));
        if numel(Pks)==3;y3=squeeze(x3(:,ps,:));end
        
        F = figure('Position', [680 648 418*1.5 330]);
        bandlims=2; 
        for nloop=1:7
            ylim([-1.4000    1.4000]);
            for pl=1:numel(Pks)
               eval(['y=y',num2str(pl),';']); 
                y_ = mean(y(61-(0:bandlims)-bandlims*(nloop-1),:),1);
                eval(['y',num2str(pl),'_=y_;']);
                if pl==1; [h1, ~] = raincloud_plot(y_(:), 2*nloop+0.5*(pl-1),'dot_dodge_amount',01,'box_on',1,'alpha',1,'bxcl',[1 0 0]);end
                if pl==2; [h1, ~] = raincloud_plot(y_(:), 2*nloop+0.5*(pl-1),'dot_dodge_amount',01,'box_on',1,'alpha',1,'bxcl',[0 0 0]);end
                if pl==3; [h1, ~] = raincloud_plot(y_(:), 2*nloop+0.5*(pl-1),'dot_dodge_amount',01,'box_on',1,'alpha',1,'bxcl',[1 1 1]/2);end
                delete(h1{1});delete(h1{2});
                hold on;
                %if pl==numel(Pks);
                %    if numel(Pks)==2; [~,p]=ttest2(y1_(:),y2_(:));end
                %    if numel(Pks)==3; [~,p]=ttest2(y1_(:),[y2_(:);y3_(:)]);end
                %    plot([1.3 1.3],2*nloop+[0 1.6]-1,'k');text(1.5,2*nloop-1,['p=',num2str(p,'%0.3f')]);
                %end
                if pl==2; [~,p]=ttest2(y1_(:),y2_(:));plot([1.0 1.0],2*nloop+[0 0.8]-1,'k');text(1.2,2*nloop-1,['p=',num2str(p,'%0.3f')]);end
                if pl==3; [~,p]=ttest2(y1_(:),y3_(:));plot([1.3 1.3],2*nloop+[0 1.6]-1,'k');text(1.5,2*nloop-1,['p=',num2str(p,'%0.3f')]);end
            end
        end
        ylim([0.05 1.6]*10);xlim([-1.5 1.5])
        view(90,-90)
        xtickformat('%0.1f')
        set(gca,'box','off')
        set(gca,'tickdir','out','ticklength',[0.03 1],'linewidth',1.5,'fontsize',15);
        
        dest = fullfile('H:\ECoG\Submit\v5\data\Sensitivity\',['Validation_',Type{bs}]);
        if ~exist(dest,'dir');mkdir(dest);end
        print(F,fullfile(dest,['Fig_',num2str(ps,'%03d'),'_',ROInames{ps},'.tiff']),'-dtiff','-r600');
        
        close all
    end
end


