clc;clear
close all
codepath = 'F:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'F:\ECoG\pca\';

load(fullfile('F:\ECoG\dcc_yyl\','states_1s.mat'));

GS_S1=[]; GS_S2=[]; 
GS_S3=[];

CutBins = 20;
CutBins_2nd = 20;

for idx=1:46
    
    
    load(fullfile('F:\ECoG\state_check\',['s',num2str(idx),'.mat']));

    %xState = ALL(:,idx);
    xState =  all_state;
    
    
    xState(xState==2)=1;
    State = interp1(1:numel(xState),xState,1.5:2:numel(xState),'nearest')';
    
    dS = State-[State(1);State(1:end-1)];
    
    load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta(:,numel(State)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(2,.3);
    beta = filtfilt(b,a,beta')';
    
    X=(beta(:,State==1)); NUM = floor(size(X,2)/100); Y=[];
    for nl=1:NUM;Y(:,nl)=mean(X(:,(nl-1)*100+(1:100)),2);end
    GS_S1 = cat(2,GS_S1,X);
    
    X=(beta(:,State==3)); NUM = floor(size(X,2)/100); Y=[];
    for nl=1:NUM;Y(:,nl)=mean(X(:,(nl-1)*100+(1:100)),2);end
    GS_S2 = cat(2,GS_S2,X);

    X=(beta(:,State==4)); NUM = floor(size(X,2)/100); Y=[];
    for nl=1:NUM;Y(:,nl)=mean(X(:,(nl-1)*100+(1:100)),2);end
    GS_S3 = cat(2,GS_S3,X);


end

StateName = {'awake';'nrem';'rem'};
ColorVec = {[238,174,81];[85,168,108];[80,98,155]};
for ss=[1 2 3]
    
    eval(['X=GS_S',num2str(ss),';']);
    Xv = squeeze(X(1:30,:,:)); 
    Ms(:,ss) = mean(Xv,2);
    STD(:,ss) = std(Xv,0,2);
    
    
end


TransitionName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};


GS_S1=[]; GS_S2=[]; 
GS_S3=[]; GS_S4=[];
St_S1=[]; St_S2=[]; 
St_S3=[]; St_S4=[];

for idx=1:46
    

    load(fullfile('F:\ECoG\state_check\',['s',num2str(idx),'.mat']));

    %xState = ALL(:,idx);
    xState =  all_state;
    
    xState(xState==2)=1;
    State = interp1(1:numel(xState),xState,1.5:2:numel(xState),'nearest')';
    Ss = State;
    
    dS = Ss-[Ss(1);Ss(1:end-1)];
    load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta(:,numel(dS)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(2,.3);
    beta = filtfilt(b,a,beta')';
    
    
    T = find(dS==+2); % state 1
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==1&unique(Ss(T(tl)+(1:CutBins_2nd)))==3
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S1 = cat(3,GS_S1,Gc(:,:,mean(St,1)~=0));
    St_S1 = cat(2,St_S1,St(:,mean(St,1)~=0));

    T = find(dS==-2); % state 2
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==3&unique(Ss(T(tl)+(1:CutBins_2nd)))==1
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S2 = cat(3,GS_S2,Gc(:,:,mean(St,1)~=0));
    St_S2 = cat(2,St_S2,St(:,mean(St,1)~=0));
    
    
    T = find(dS==+1); % state 3
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==3&unique(Ss(T(tl)+(1:CutBins_2nd)))==4
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S3 = cat(3,GS_S3,Gc(:,:,mean(St,1)~=0));
    St_S3 = cat(2,St_S3,St(:,mean(St,1)~=0));

    
    T = find(dS==-3); % state 4
    T (T<=CutBins | T >=numel(Ss)-CutBins) = [];
    Gc = zeros(100,CutBins*2+1,numel(T));
    St = zeros(CutBins*2+1,numel(T));
    for tl=1:numel(T)
        if unique(Ss(T(tl)+(-CutBins_2nd:-1)))==4&unique(Ss(T(tl)+(1:CutBins_2nd)))==1
            Gc(:,:,tl) = beta(:,T(tl)+(-CutBins:CutBins));
            St(:,tl) = Ss(T(tl)+(-CutBins:CutBins));
        end
    end
    GS_S4 = cat(3,GS_S4,Gc(:,:,mean(St,1)~=0));
    St_S4 = cat(2,St_S4,St(:,mean(St,1)~=0));

    
end



TransitionName = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};
PC_ColorVec = {[254,154,118];[222,90,0];[103,154,155];[99,115,115]};

F = figure('Position', [680 572 997 406]);
for ss=1:numel(TransitionName)
    
    x_ = (-CutBins:CutBins)*2;

    if ss==1;Bs=[1 2];end
    if ss==2;Bs=[2 1];end
    if ss==3;Bs=[2 3];end
    if ss==4;Bs=[3 1];end
    
    eval(['GS=GS_S',num2str(ss),';']);
    eval(['St=St_S',num2str(ss),';']);
    
    clear Dis*
%     Gx = mean(GS(1:30,:,:),3);
    Gx = GS(1:30,:,:);
    bs = 100;
    if ss==3 | ss==4; bs=20;end
    
    for iz=1:floor(size(Gx,3)/bs)
        for iy = 1:size(Gx,2)
            Dis1(iy,iz) = sqrt(sum((Ms(:,Bs(1))'-mean(Gx(:,iy,iz:floor(size(Gx,3)/bs):end),3)').^2));
            Dis2(iy,iz) = sqrt(sum((Ms(:,Bs(2))'-mean(Gx(:,iy,iz:floor(size(Gx,3)/bs):end),3)').^2));
        end
    end
    y = mean(Dis1,2)-mean(Dis2,2);
    
    vars{1} = 1:numel(y);
    pars(1) = 30;        lb(1) = pars(1)-200;    ub(1) = pars(1)+200;
    pars(2) = .21;       lb(2) = pars(2)*0.01;   ub(2) = pars(2)*100;
    pars(3) = 0.1;       lb(3) = 0;             ub(3) = pars(3)*100;
    pars(4) = 1;       lb(4) = 0;             ub(4) = pars(4)*100;
        
    options = optimset('LargeScale','on',...
        'Algorithm', 'trust-region-reflective',...
        'TolFun',10^(-100),'TolX',10^(-10),...
        'MaxFunEvals',100000,'MaxIter',100000);
        
    func = @(pars,vars) pars(3)*gamcdf(vars{1},pars(1),pars(2))-pars(4);
    [pars0,~,~,~] = lsqcurvefit(func,pars,vars,y(:)',lb,ub,options);
    fitted_BOLD = pars0(3)*gamcdf(vars{1},pars0(1),pars0(2))-pars0(4);

    
    subplot(2,4,ss);
    plot(x_',Dis1','o','MarkerSize',4,'MarkerFaceColor',ColorVec{Bs(1)}/255,'MarkerEdgeColor','none'); hold on;
    plot(x_',Dis2','o','MarkerSize',4,'MarkerFaceColor',ColorVec{Bs(2)}/255,'MarkerEdgeColor','none'); hold on;
    %legend({['Dis to ',StateName{Bs(1)}];['Dis to ',StateName{Bs(2)}]},...
    %    'box','off','location','best');
    xlim([min(x_),max(x_)])
    title([replace(TransitionName{ss},'_',' 2 ')])
    set(gca,'box','off','linewidth',1,'fontsize',12,'tickdir','out','ticklength',[.2 .1]/4);
    
    subplot(2,4,ss+4);
    plot(x_,y);hold on;
    plot(x_,fitted_BOLD,'k');
    xlim([min(x_),max(x_)])

    plot(x_,gradient(fitted_BOLD));
    TTP = x_(gradient(fitted_BOLD)==max(gradient(fitted_BOLD)))
    [~,LOCS] = findpeaks( -abs(gradient(fitted_BOLD)-max(gradient(fitted_BOLD))/2));
     FWHM = x_(LOCS(2)) - x_(LOCS(1))
    set(gca,'box','off','linewidth',1,'fontsize',12,'tickdir','out','ticklength',[.2 .1]/4);
    1;
end