clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\Activity_transition\';

State = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};

%% global signal transition
filepath = 'H:\ECoG\G_S\G_S\';
load('H:\ECoG\dcc_yyl\states_tr.mat');
GS_S1=[]; GS_S2=[]; 
GS_S3=[]; GS_S4=[]; 

CutBins = 30;

for idx=1:46
    cd(fullfile(filepath,num2str(idx,'%02d')));
    GS = load('G_S.txt');
    Ss = ALL(:,idx);
    Ss(Ss==2)=1;
    dS = Ss-[Ss(1);Ss(1:end-1)];
    
    T = find(dS==+2); % state 1
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T);Gc(:,tl) = GS(T(tl)+(-CutBins:CutBins));end
    GS_S1 = cat(2,GS_S1,Gc);
    
    T = find(dS==-2); % state 2
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T);Gc(:,tl) = GS(T(tl)+(-CutBins:CutBins));end
    GS_S2 = cat(2,GS_S2,Gc);
    
    
    T = find(dS==+1); % state 3
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T);Gc(:,tl) = GS(T(tl)+(-CutBins:CutBins));end
    GS_S3 = cat(2,GS_S3,Gc);

    
    T = find(dS==-3); % state 4
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T);Gc(:,tl) = GS(T(tl)+(-CutBins:CutBins));end
    GS_S4 = cat(2,GS_S4,Gc);

    1;
end

shift = -CutBins:CutBins;
F = figure;
subplot(2,2,1);
plot(shift*2,mean(GS_S1,2)); 
ylim([-.1 .1]); xlim([-60 60]);
hold on;plot([0  0],[-.1 .1]*4);plot([-60 60],[0 0]);
subplot(2,2,2);
plot(shift*2,mean(GS_S2,2));
ylim([-.1 .1]); xlim([-60 60]);
hold on;plot([0  0],[-.1 .1]*4);plot([-60 60],[0 0]);
subplot(2,2,3);
plot(shift*2,mean(GS_S3,2));
ylim([-.1 .1]*4); xlim([-60 60]);
hold on;plot([0  0],[-.1 .1]*4);plot([-60 60],[0 0]);
subplot(2,2,4);
plot(shift*2,mean(GS_S4,2));
ylim([-.1 .1]*4); xlim([-60 60]);
hold on;plot([0  0],[-.1 .1]*4);plot([-60 60],[0 0]);
%% global relative time lag map
%

close all
for sl=1:numel(State)
    load(fullfile(WholePath,[State{sl},'_modified.mat']));
    
    %y = load(fullfile(WholePath,[State{sl},'_x.mat']));
    %ShiftTime = -y.ans;
    
    %
    shift = -CutBins:CutBins;
%     Labs = {'Isocortex_primary';'Isocortex_HigherOrder';'HPF';'STR';'PAL';'TH';'HY';'MB'};
%     Defs = {[1:18];[18:31];[32:36];[38:41];[42:45];[46:55];[56:59];[60:69]};
%     CV = {[39,143,87];[39,143,87]*.7;[123,75,30];[238,173,77];[183,28,37];[220,0,0];[153,0,255];[154,27,91]};
    Labs = {'VIS';'AUD';'MO';'SS';'ACA';'RSC';'PFC';'HPF';'STR';'PAL';'TH';'HY';'MB'};
    Defs = {[15:18];[12:14];[1:2];[3:9];[19:20];[25:27];[21:24];[32:36];[38:41];[42:45];[46:55];[56:59];[60:69]};
    %[~,s] = sort(X);
    Ymatrix = test;%YY(s,:);
    
    filePATH = fullfile(WholePath,'tran_61',State{sl});
    for lp=1:61
        hdr = spm_vol(fullfile(filePATH,['avg',num2str(lp),'.nii']));
        img = spm_read_vols(hdr);
        if lp==1;fMRI_4D=zeros([size(img),61]);end
        fMRI_4D(:,:,:,lp)=img;
    end
    
    filepath = 'H:\ECoG\dcc_yyl\';
    Labels_NII = fullfile(filepath,'Label_69.nii');
    Labels_Excel = fullfile(filepath,'Label_69.xlsx');
    [~,~,CellData] = xlsread(Labels_Excel);
    ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
    MajorRegion = cat(1,ExpTable.Major_Region);
    SubRegion = cat(1,ExpTable.Sub_Region);
    MinRegion = cat(1,ExpTable.Min_Region);
    L = spm_read_vols(spm_vol(Labels_NII));

    MatSer = zeros(numel(Labs),size(test,2));
    for tl=1:numel(Labs)
        Mask = L>=min(Defs{tl}) & L<=max(Defs{tl});
        V = fmask(fMRI_4D,Mask);
        mC = nanmean(V,1);
        [fa,fb]=butter(3,.3);
        %mC = filtfilt(fa,fb,mC);
        MatSer(tl,:) = mC;
    end
    
    L1 = cell(numel(label1),1);
    for xl=1:numel(label1)
        pin = cellfun(@(x) strcmp(x,label1{xl}), ExpTable.abbre);
        L1(xl) = SubRegion(pin);
    end
    
    List(:,3) = label1(:);
    List(:,2) = repmat({1},numel(1:69),1);
    List(:,1) = L1;
%     SubUnique = unique(MinRegion);%setdiff(unique(MinRegion),unique(SubRegion));
%     for sul = 1:numel(SubUnique)
%         if ~contains(SubUnique(sul),'Others')
%             pin = cellfun(@(x) strcmp(x,SubUnique(sul)),MinRegion);
%             NUMB = numel(find(pin==1));
%             ParentRegion = unique(SubRegion(pin));
%             ListPlus = [ParentRegion,{NUMB},SubUnique(sul)]; 
%             List = cat(1,ListPlus,List);
%         end
%     end
    
    Ls = (List);
    colorList = jet(69);
    F = figure('Position', [680 443 548 535]);
    subplot(1,2,1)
    sankeyHdl=sankey2([],'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'List',Ls,'Color',colorList);
    axis off
    subplot(1,2,2)
    gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
    drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
    defaultMap = [gdmap;flipud(drmap);];
    %YY(YY==0)=nan;
    imagesc(shift*2,1:size(test,1),test);colormap(defaultMap);
    caxis([-.2 .2])
    if sl==3|sl==4;caxis([-.6 .6]);end
   % saveas(F,fullfile(WholePath,[State{sl},'_series.emf']));

%     F = figure;
%     subplot(5,1,1);
%     eval(['GS_=mean(GS_S',num2str(sl),',2);']);
%     plot(shift*2,GS_); ylim([-.1 .1]*2); xlim([-60 60]);
%     if sl>=3;ylim([-.4,.4]);end
%     set(gca,'box','off');
%     subplot(5,1,2:5);
%     plot(MatSer'+(1:numel(Labs))/10);
end    
    
    %}

    
 %{   
% % % % % %     clear CC GS Sd
% % % % % %     eval(['GS_=mean(GS_S',num2str(sl),',2);']);
% % % % % %     vars = GS_(:)';
% % % % % %     %vars = [zeros(1,30),.5,ones(1,30)];
% % % % % %     
% % % % % %     
% % % % % %     shift = -CutBins:CutBins;
% % % % % %     for shl=1:numel(shift)
% % % % % %        CC(:,shl) = corr(YY',circshift(vars(:),shift(shl),1)); 
% % % % % %     end
% % % % % %     
% % % % % %     
% % % % % %     for cl=1:size(CC,1)
% % % % % %         x = CC(cl,:); y = find(x==max(x));
% % % % % %         Sd(cl)=y;
% % % % % %         ShiftTime(cl) = shift(y)*2;
% % % % % %     end
% % % % % %     [~,s] = sort(Sd);
%
    X0 = [X;ShiftTime]';
    Tim = sortrows(X0,'ascend');

    %figure;imagesc(YY(s,:));
    
    filepath = 'H:\ECoG\dcc_yyl\';
    Template_NII = fullfile(filepath,'Template_Mouse_X20.nii');
    Labels_NII = fullfile(filepath,'wLabel_69.nii');
    Labels_Excel = fullfile(filepath,'Label_69.xlsx');

    
    [~,~,CellData] = xlsread(Labels_Excel);
    ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
    MajorRegion = cat(1,ExpTable.Major_Region);
    ihdr = spm_vol(Template_NII);
    T = spm_read_vols(spm_vol(Template_NII));
    L = spm_read_vols(spm_vol(Labels_NII));
    
    y = Tim(:,2);
%     y(y>=0) = y(y>=0)+4;
%     y(y<0) = y(y<0)-4;
    
    I = zeros(size(T));
    for ilx=1:numel(y)
        J = zeros(size(T));
        J(L==ilx) = y(ilx);
        I = I+J;
    end
    I = I+flip(I,1);
    I = smooth3(I,'box',11);
    I(L==0)=0;
    hdr = ihdr;  hdr.dt=[16,0];
    hdr.fname = fullfile(WholePath,['X_',State{sl},'_timelag.nii']);
    spm_write_vol(hdr,I);
    
    Func3D = fullfile(WholePath,['X_',State{sl},'_timelag.nii']);
    delete(fullfile(WholePath,['n','X_',State{sl},'_timelag.nii']));
    [f1_rgb,f2_rgb,f3_rgb]=Colormap_3Dviewer_v2(Func3D,Template_NII,[-30 -.1 .1 30]/2,codepath,1);

    I = cat(1,f2_rgb.cdata,f3_rgb.cdata,f1_rgb.cdata);
    for i=1:3; I(:,:,i) = medfilt2(I(:,:,i),[3 3]);end
    J = rgb2hsv(I);J(:,:,3)=J(:,:,3)*1.2;
    K = hsv2rgb(J);
    imwrite(uint8(K),fullfile(WholePath,[State{sl},'_3D.tiff']),'Resolution',300);
    
   
     close all       
    1;

end
%}
%{
close all

for sl=1:numel(State)
    load(fullfile(WholePath,[State{sl},'_raw.mat']));
    
    if sl==1;YY(:,[1:10,end-9:end])=[];end
    vars{1} = 1:size(YY,2);
    figure('position',[10 10 1600 900]);
    for idx=1:size(YY,1)
        CC = corr(YY(idx,:)',vars{1}');
        Y_ = YY(idx,:) * CC/abs(CC);
        [bx,ax] = butter(3,.2);
        actual_BOLD = filtfilt(bx,ax,Y_);
        pars(1) = 30;        lb(1) = pars(1)-200;    ub(1) = pars(1)+200;
        pars(2) = .21;       lb(2) = pars(2)*0.01;   ub(2) = pars(2)*100;
        pars(3) = 0.1;       lb(3) = 0;             ub(3) = pars(3)*10;
        pars(4) = 0.1;       lb(4) = 0;             ub(4) = pars(4)*10;
        
        options = optimset('LargeScale','on',...
            'Algorithm', 'trust-region-reflective',...
            'TolFun',10^(-100),'TolX',10^(-10),...
            'MaxFunEvals',100000,'MaxIter',100000);
        
        
        
        func = @(pars,vars) pars(3)*gamcdf(vars{1},pars(1),pars(2))-pars(4);
        [pars0,~,~,~] = lsqcurvefit(func,pars,vars,actual_BOLD,lb,ub,options);
        fitted_BOLD = pars0(3)*gamcdf(vars{1},pars0(1),pars0(2))-pars0(4);
        
        subplot(7,10,idx);
        plot(YY(idx,:)); 
        if sl==1 | sl==2 ; ylim([-0.1 .1]*2);end
        if sl==3 | sl==4 ; ylim([-0.1 .1]*6);end
        title(label1{idx});
        hold on;plot(fitted_BOLD* CC/abs(CC),'linewidth',2);
        stats = regstats(fitted_BOLD'*CC/abs(CC),Y_','linear');
        p(idx) = stats.tstat.pval(2);
        R(idx) = stats.rsquare;
    end
end
%}