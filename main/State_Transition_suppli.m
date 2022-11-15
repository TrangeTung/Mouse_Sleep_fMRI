clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'H:\ECoG\';

State = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};

%% State transition motif
filepath = 'H:\ECoG\Activity_transition\tran_31\';
Temp3D = fullfile('H:\ECoG\ripple-average\','Template_Mouse_v38.nii');
% for sl=1:numel(State)
% Nmean = 1;I=[];
% barlim = [-.1 -.02 .02 .1]*5;
% if sl==1;barlim = [-.1 -.02 .02 .1]*2;end
% if sl==2;barlim = [-.1 -.02 .02 .1]*2;end
% for f = 2:2:31;
%     Func3D = fullfile(filepath,State{sl},['avg',num2str(f),'.nii']);
%  Img_RGB=Colormap(   'statfile',Func3D,...
%                     'bgmfile',Temp3D,...
%                     'slice',9:2:27,...
%                     'bar_value',barlim,...
%                     'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
%                     'dest',fullfile(filepath,State{sl}),...
%                     'mapname',['s'],...
%                     'cluster',10);
%       
%     I = cat(1,I,Img_RGB);
% 
% 
% 
% end
% imwrite(uint8(I),fullfile(filepath,State{sl},'Transition.tiff'));
% end


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



v = VideoWriter(fullfile(WholePath,'Activity_transition','StateTransition_1.avi'));
v.FrameRate=1;
open(v); 

close all
for f = 1:61

    FrameI = cell(1,4);
for sl=1:numel(State)
    
    
        
    filePATH = fullfile(WholePath,'Activity_transition','tran_61',State{sl});
    Func3D = fullfile(filePATH,'Transition.nii');
    
    Temp3D = fullfile(WholePath,'2nd','Template_Mouse_v38.nii');
    barlim = [-.101 .101];
    if sl==3|sl==4;barlim = [-.301 .301];end
    [f1_rgb,f2_rgb,f3_rgb,f4_rgb]=Colormap_3Dviewer_v1(Func3D,Temp3D,barlim,codepath,f);
    
    Ix = f2_rgb.cdata(81:350,71:490,:);
    Iz = f3_rgb.cdata(:,:,:);
    Iy = f1_rgb.cdata(:,:,:);
    Iwx = f4_rgb.cdata;%(:,41:540,:);
    Iw = imresize(Iwx,size(Iz,1)/size(Iwx,1)*.9);
    
    I2 = cat(1,Iz);
    I1 = imresize(Iwx,size(I2,2)/size(Iwx,2));
    I0 = cat(1,I1,I2);
    I3 = 255*ones(size(I0));
    I3(200+(1:size(Iy,1)),:,:)=Iy;
    
    %I1 = cat(2,Iy,permute(Ix,[2 1 3]));
    %I2 = 255*ones(size(Iz,1),size(I1,2),3);
    %I2(:,1:size(Iz,2),:)=Iz;
    %I2(end-size(Iw,1)+1:end,(end+1-size(Iw,2)):end,:)=Iw;
    
    I = cat(2,I0,I3);
    
    %{
    if f==0
        load(fullfile(WholePath,'Activity_transition',[State{sl},'_raw.mat']));
        y = load(fullfile(WholePath,'Activity_transition',[State{sl},'_x.mat']));
        ShiftTime = -y.ans;
        
        CutBins = 30;
        shift = -CutBins:CutBins;
        
        F = figure('Position', [1349 213 320 751],'color',[1 1 1]);
        imagesc(shift*2,1:size(YY,1),YY);
        set(gca,'Ytick',1:size(YY,1),'Yticklabel',label1);
        xlabel('Time (s)');axis tight
        set(gca,'fontsize',10);
        caxis([-.201 .201]);colormap('jet');
        colorbar('NorthOutside');
        if sl==3|sl==4;caxis([-.601 .601]);end
            
        Q = getframe(F);
        Qc = Q.cdata;
        Qcr = imresize(Qc,size(I,1)/size(Qc,1));
        
        if sl==1;QcrT{sl,1}=insertText(Qcr,[0,0],'AW -> NREM','fontsize',30);end
        if sl==2;QcrT{sl,1}=insertText(Qcr,[0,0],'NREM -> AW','fontsize',30);end
        if sl==3;QcrT{sl,1}=insertText(Qcr,[0,0],'NREM -> REM','fontsize',30);end
        if sl==4;QcrT{sl,1}=insertText(Qcr,[0,0],'REM -> AW','fontsize',30);end
    end
    %}
    
    
    
    %if sl==1;FrameI{1,2}=cat(2,I,QcrT{sl,1});end
    %if sl==2;FrameI{2,2}=cat(2,I,QcrT{sl,1});end
    %if sl==3;FrameI{1,1}=cat(2,I,QcrT{sl,1});end
    %if sl==4;FrameI{2,1}=cat(2,I,QcrT{sl,1});end
    
    
    
    
    
    CutBins = 30;
    shift = -CutBins:CutBins;
    F = figure('Position',  [680 397 598 581],'color',[1 1 1]);
    eval(['GS=mean(GS_S',num2str(sl),',2)*100;']);
    subplot('position',[.25 .37 .6 .6])
    plot(shift*2,GS,'linewidth',2,'color',[.7 .7 .7]);
    xlabel('Time (s)'); ylabel('Global signal (%)');
    hold on;
    plot(shift(1:f)*2,GS(1:f),'linewidth',6,'color',[1.0 .0 .0]);
    plot([0 0],[-100 100],'k--');
    plot([-60 60],[0 0],'k--');
    set(gca,'box','off','linewidth',1.5,'tickdir','out','fontsize',15);
    ylim([-3 3]);
    if sl==3 | sl==4;ylim([-30 30]);end
    subplot('position',[.25 .05 .6 .3]);
    axis off
    colormap('jet'); caxis(barlim*10.0001);
    c=colorbar('South','FontSize',18,'linewidth',1.5);
    c.Label.String = 'BOLD signal (z)';
    
    Q = getframe(F);
    Qc = Q.cdata;
    Qcr = imresize(Qc,size(I,2)/size(Qc,2),'bilinear');
    Qcrm = cat(1,255*ones(100,size(Qcr,2),3),Qcr);
    close all
    
    if sl==1;J=insertText(Qcrm,[0,0],'AW -> NREM','fontsize',30);end
    if sl==2;J=insertText(Qcrm,[0,0],'NREM -> AW','fontsize',30);end
    if sl==3;J=insertText(Qcrm,[0,0],'NREM -> REM','fontsize',30);end
    if sl==4;J=insertText(Qcrm,[0,0],'REM -> AW','fontsize',30);end
    J = insertText(J,[0,size(J,2)*.9],['Time: ',num2str((f-31)*2),' s'],...
        'fontsize',30,'BoxColor','white');
    
    FrameI{1,sl}=cat(1,J,I);
    
    1;
    
    %{
    filePATH = fullfile(WholePath,'tran_61',State{sl});
    Iy = [];
    for lp=31-6:1:31+6
        filename = fullfile(filePATH,['avg',num2str(lp),'.nii']);
        Func_Img_3D = spm_read_vols_4D(spm_vol(filename));
        Func_Img_3D = Func_Img_3D/2+flip(Func_Img_3D,1)/1.5;
        Func_Img_3D = smooth3(Func_Img_3D,'box',3);
        slice = 13:2:25;
        Ix = [];
        for sl=1:numel(slice)
            I = MY_display_function_map_3D_nothreshold(Func_Img_3D,[-.2 .2],slice(sl));
            Ix = cat(1,Ix,I);
        end
        
        Iy = cat(2,Iy,Ix);
    end
    %}
    
    
end

W = cell2mat(FrameI);
[x,y,z] = size(W);
W(:,round(y/4*1)+(-1:1),:) = 100;
W(:,round(y/4*2)+(-1:1),:) = 100;
W(:,round(y/4*3)+(-1:1),:) = 100;

Wr=W;
%Wr = imresize(W,[1080 1920],'bilinear');

writeVideo(v,Wr);

end
close(v);




for speed = 1:3

Type = {'awake_nrem';'nrem_awake';'nrem_rem';'rem_awake'};
Name = {'AW to NREM';'NREM to AW';'NREM to REM';'REM to AW'};

v = VideoWriter(fullfile(WholePath,['Transition_slice_X',num2str(speed),'speed.avi']));
v.FrameRate=speed;
open(v);


for f=1:61

    I = [];

    for tl = 1:numel(Type)

        filePATH = fullfile(WholePath,'Activity_transition','tran_61',Type{tl});
        Func3D = fullfile(filePATH,'Transition.nii');


        close all
        Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');

        Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                        'bgmfile',Temp3D,...
                        'slice',9:2:27,...
                        'bar_value',[-.1 -.02 .02 .1]*3,...
                        'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                        'dest',fullfile(WholePath),...
                        'mapname',['s'],...
                        'cluster',10);
        Ix = uint8(Img_RGB);
        Iy = uint8(zeros([130 1140 3]));
        %if tl==1;
        Iy = insertText(Iy,[size(Ix,2)*0.4,0],['Time: ',num2str((f-31)*2),' s'],...
            'fontsize',25,'BoxColor','black','TextColor','yellow');
        %end
        Iy = insertText(Iy,[size(Ix,2)*0,0],Name{tl},...
            'fontsize',30,'BoxColor','black','TextColor','white');

        Iy(end-size(Ix,1)+1:end,:,:) = Ix;
        
        I = cat(1,I,Iy);
        
    end
    
    writeVideo(v,I);

    
    
    
end
    close(v);



    end
