clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\pca\';

%% Group PCA
%{
COEFF = [];
for idx=1:46
    load(fullfile(WholePath,num2str(idx,'%02d'),'pca_t.mat'));
    COEFF = cat(2,COEFF,coeff);
end
[Gcoeff, Gscore, ~, ~, explained, ~] = pca(COEFF','NumComponents',100);
save(fullfile(WholePath,'GroupWeight.mat'),'Gcoeff','Gscore','-v7.3');

load(fullfile(WholePath,'GroupWeight.mat'));
for idx=30:46
    
    load(fullfile(WholePath,num2str(idx,'%02d'),'pca_t.mat'));
 
    
    %
    tic;
    Gs = coeff*score';
    
    V = Gs; clear Gs
    DM = Gcoeff;
    lamda = 10:-1:-20;
    fa = 10^lamda(13);
    beta = pinv(DM'*DM+fa*eye(size(DM,2)))*(DM'*V);

    Vec = (V(:)');
    R2_gpu = (zeros(size(beta,1),1));
    %p = parpool(4);
    for sloop = 1:size(beta,1)
        VecP0 = (beta(sloop,:)'*DM(:,sloop)')';
        VecP = VecP0(:)';  VecP0=[];
        mu = mean(VecP,2).*mean(Vec,2)-mean(Vec.*VecP,2);
        md = mean(VecP,2).^2 - mean(VecP.^2,2);
        b = mean(Vec,2) - mu./md.*mean(VecP,2);
        Vecf = mu./md.*VecP+b;               VecP =[];
        SSR = sum((Vecf-mean(Vec,2)).^2,2);  Vecf =[];
        SST = sum((Vec-mean(Vec,2)).^2,2);
        R2_gpu(sloop,1) = SSR./SST;
    end
    %delete(p);
    R2 = (R2_gpu);
    clear Vec
                
    save(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'),'beta','R2','-v7.3');
    toc;
    %
end
%}


%% Quantitative gPCA
%{
load(fullfile(WholePath,'GroupWeight.mat'));
ihdr = spm_vol(fullfile(WholePath,'mask_RAM.nii'));
lmask = spm_read_vols(ihdr);
fMRI_4D = funmask(Gcoeff(:,1:10),lmask)*100;
for loop=1:10
   Xhdr(loop)=ihdr;Xhdr(loop).n=[loop,1];
   Xhdr(loop).dt=[16,0];
   Xhdr(loop).fname = fullfile(WholePath,'PCweight.nii');
end
spm_write_vol_4D(Xhdr,fMRI_4D);



%}

%{
Func3D = fullfile(WholePath,'PCweight.nii');

Temp3D = fullfile('H:\ECoG\','Activation','2nd','mean_T2.nii');
barlim = [-.601 .501];
for f=1:10
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
%     
%     Ix = f2_rgb.cdata(81:350,71:490,:);
%     Iz = f3_rgb.cdata(81:350,101:520,:);
%     Iy = f1_rgb.cdata(:,101:460,:);
%     Iwx = f4_rgb.cdata(:,41:540,:);
%     Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);
% 
%     I1 = cat(2,Iy,permute(Ix,[2 1 3]));
%     I2 = 255*ones(size(Iz,1),size(I1,2),3);
%     I2(:,1:size(Iz,2),:)=Iz;
%     I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
%     
%     I = cat(1,I1,I2);
%     
%     I1 = double(f4_rgb.cdata);
%     I2 = double(f3_rgb.cdata);
%     
%     I = zeros(720,760,3);
%     I1(end-100+1:end,end-100+1:end,:)=0;
%     I(1:size(I1,1),1:size(I1,2),:) = I(1:size(I1,1),1:size(I1,2),:) + double(I1);
%     I2(1:100,1:250,:)=0;
%     I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:) = I(end-size(I2,1)+1:end,end-size(I2,2)+1:end,:)+double(I2);
%     I(1:end-size(I2,1)+1,size(I1,2)+1:end,:) = 255;
%     I(size(I1,1)+1:end,1:end-size(I2,2)+1,:) = 255;

    I1 = double(f1_rgb.cdata);
    I2 = double(f4_rgb.cdata);
    I3 = double(f3_rgb.cdata);
    I = cat(1,cat(1,I1,I2),I3);
    dest = fullfile(WholePath,['PC_',num2str(f,'%02d')]);
    mkdir(dest);
    
    imwrite(uint8(I),fullfile(dest,['PC_',num2str(f,'%02d'),'.tiff']));
    close all;
    
%     Colormap(   'statfile',[Func3D,',',num2str(f)],...
%                 'bgmfile',Temp3D,...
%                 'slice',9:2:27,...
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

%% Individual  trajectory
%{
load(fullfile('H:\ECoG\dcc_yyl\','states_tr.mat'));
X=[];
F = figure;
for idx=1:46
    load(fullfile('H:\ECoG\dcc_yyl',num2str(idx,'%02d'),'timeseries_138.mat'));%GroupPCs
    X = cat(2,X,roi);
end

% ref
LS = spm_read_vols(spm_vol('H:\ECoG\dcc_yyl\Label_69.nii'));
PCw = spm_read_vols(spm_vol('H:\ECoG\pca\PCweight.nii'));
Vref = zeros(138,10);
for ix=1:138
   Vref(ix,:) = mean(fmask(PCw,LS==ix),1);
end



[coeff,score,latent, tsquared,explained] = pca([X'],'NumComponents',20);%,'Algorithm','als','Coeff0' ,Vref


y(:,1) = explained(1:100)/sum(explained(1:100))*100;
y(:,2) = cumsum(y(:,1));

dest = fullfile('H:\ECoG\','GroupPCA');
%save(fullfile(dest,'GroupPCA.mat'),'coeff','score');



%% Exp R2





%}

load(fullfile('H:\ECoG\dcc_yyl\','states_tr.mat'));
ns = 300;xlc = 0;

Sy=[];St=[];
for idx=1:46
    State = ALL(:,idx);
    State(State==2)=1;
    uniqueS = unique(State);
    
    dS = State-[State(1);State(1:end-1)];
    
    load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta(:,numel(State)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(3,.1);
    beta = filtfilt(b,a,beta')';
    P = double(dS~=0);
    
    
    [~,LOCS] = findpeaks(P,'MinPeakHeight',0.5,'MinPeakDistance',10);
    
    
    for xl=1:numel(LOCS)-1
        
        xlc=xlc+1;
        S = beta(1:10,LOCS(xl):LOCS(xl+1));
        Si = interp1((1:size(S,2))/size(S,2)*ns, S', 1:ns, 'nearest')';
        Si = fillmissing(Si','nearest')';
        Sy = cat(2,Sy,Si);
        
        
        S = State(LOCS(xl):LOCS(xl+1))';
        Si = interp1((1:size(S,2))/size(S,2)*ns, S', 1:ns, 'nearest')';
        Si = fillmissing(Si','nearest')';
        St = cat(2,St,Si(:)');
        
        
    end
end




Sall = [];%zeros(10,ns,3600);
Stall = [];%zeros(10,ns,3600);
Stateall = [];%zeros(ns,3600);
Transall = [];%zeros(ns,3600);
S_MouseIdx = [];%zeros(3600,1);
T_MouseIdx = [];%zeros(3600,1);
xlc=0;xlt=0;

dS = St-[St(1),St(1:end-1)];
P = double(dS~=0);
[~,LOCS] = findpeaks(P,'MinPeakHeight',0.5,'MinPeakDistance',10);

for xl=1:numel(LOCS)-1
        
    xlc=xlc+1;
    S = Sy(1:10,LOCS(xl)+(51:250));
    Sall(:,:,xlc)=S;

    S = St(LOCS(xl)+(51:250))';
    Stateall(:,xlc)=S;

end


for xl=1:numel(LOCS)

    if LOCS(xl)>20 & LOCS(xl)<numel(dS)-20
        xlt=xlt+1;
        S = Sy(1:10,LOCS(xl)+(-100:100));
        Stall(:,:,xlt)=S;

        S = St(LOCS(xl)+(-100:100))';
        Transall(:,xlt)=S;

    end
end



%{
for idx=1:46
    State = ALL(:,idx);
    State(State==2)=1;
    uniqueS = unique(State);
    
    dS = State-[State(1);State(1:end-1)];
    
    load(fullfile(WholePath,num2str(idx,'%02d'),'GroupPCs.mat'));%
    beta(:,numel(State)+1:end)=[];
    beta = (beta-mean(beta,2) )./std(beta,0,2);
    [b,a] = butter(2,.2);
    beta = filtfilt(b,a,beta')';
    y=hilbert(beta(1,:));
    P=atan2(imag(y),real(y));
    P = double(dS~=0);
    
    
    [~,LOCS] = findpeaks(P,'MinPeakHeight',0.5,'MinPeakDistance',10);
    
    Sx=[];Sy=[];
    for xl=1:3%numel(LOCS)-1
        
        xlc=xlc+1;
        S = beta(1:10,LOCS(xl):LOCS(xl+1));
        Si = interp1((1:size(S,2))/size(S,2)*ns, S', 1:ns, 'nearest')';
        Si = fillmissing(Si','nearest')';
        
        Sx = cat(2,Sx,S);
        Sy = cat(2,Sy,Si);
        
%         Sall(:,:,xlc)=Si;
%         
%         S = State(LOCS(xl):LOCS(xl+1))';
%         Si = interp1((1:size(S,2))/size(S,2)*ns, S', 1:ns, 'nearest')';
%         Si = fillmissing(Si,'nearest');
%         Stateall(:,xlc)=Si;
%         
%         
%         S_MouseIdx(xlc)=idx;
    end
    
    
    for xl=1:numel(LOCS)
        
        if LOCS(xl)>20 & LOCS(xl)<numel(State)-20
            xlt=xlt+1;
            S = beta(1:10,LOCS(xl)+(-20:20));
            Si = interp1((1:size(S,2))/size(S,2)*ns, S', 1:ns, 'nearest')';
            Si = fillmissing(Si','nearest')';
            Stall(:,:,xlt)=Si;

            S = State(LOCS(xl)+(-20:20))';
            Si = interp1((1:size(S,2))/size(S,2)*ns, S', 1:ns, 'nearest')';
            Si = fillmissing(Si,'nearest');
            Transall(:,xlt)=Si;


            T_MouseIdx(xlt)=idx;
        end
    end
    
    1;
    
    
    
    
end
%}

%% Pure state attractor

F = figure;
State = {'awake';'nrem';'rem'};
for ss=[1 2 3]
   
    %X = std(Stateall,0,1)==0;
    Ys = nanmedian(Stateall,1);
    if ss==1; Y=Ys==1; ColorVec = [1,0,0];end
    if ss==2; Y=Ys==3; ColorVec = [0,1,0];end
    if ss==3; Y=Ys==4; ColorVec = [0,0,1];end
    
    %ColorVec = hsv(ns);
    Binary =  find(Y==1);
%     for idx = 1:46
%        pin  =  find(( Binary(:) & MouseIdx==idx )==1);
%        if ~isempty(pin)
%            clear S
%            S = median(Sall(:,:,pin),3);
% 
%            scatter3(S(2,:)',S(3,:)',S(4,:)','MarkerFacealpha',0.1,...
%            'MarkerEdgeColor','none','MarkerFaceColor',ColorVec);
%            hold on;
%        else
%            idx
%        end
%     end
    S = nanmean(Sall(:,:,Binary),3);
    eval(['S_',num2str(State{ss}),'=S;']);
    xyz = S(1:3,:);
    %plot3(S(1,:),S(2,:),S(3,:),'color',ColorVec,'linewidth',3);hold on;
    fnplt(cscvn(squeeze(xyz(:,[1:end ]))),3);hold on;
    
    
    St = Sall(:,:,Binary);
    dest = fullfile(WholePath,'Trajectory');
    mkdir(dest);
    save(fullfile(dest,['Trajectory_',State{ss},'.mat']),'St')
    
    
end


%% transition attractor
State = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};

for ss=1:4
    if ss==1; pin=+2; end % state 1
    if ss==2; pin=-2; end % state 2
    if ss==3; pin=+1; end % state 3
    if ss==4; pin=-3; end % state 4
    
    dS = Transall - [Transall(1,:);Transall(1:end-1,:)];
    [sx,sy]=size(dS);
    Iy = find(dS(round(sx/2),:)==pin);
    Binary = Iy;
    
    
    S = nanmean(Stall(:,:,Binary),3);
    RM = strsplit(State{ss},'_');
    eval(['Sa=S_',RM{1},';']); 
    eval(['Sb=S_',RM{2},';']); 
    Sy_st=mean(Sa(:,end-19:end),2);
    Sy_ed=mean(Sb(:,1:20),2); 
    Sx_st = mean(S(:,1:20),2);
    Sx_ed = mean(S(:,end-19:end),2);
    
    k = (Sy_ed-Sy_st)./(Sx_ed-Sx_st);
    b = Sy_ed-Sx_ed.*k;
    
    S_ = S.*k+b;
    xyz = S_(1:3,:);
    fnplt(cscvn(squeeze(xyz(:,[1:end ]))),'k',3);hold on;
    
    
    St = Stall(:,:,Binary).*k+b;
    dest = fullfile(WholePath,'Trajectory');
    mkdir(dest);
    save(fullfile(dest,['Trajectory_',State{ss},'.mat']),'St')

    
end
