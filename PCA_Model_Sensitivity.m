clc;clear
close all
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'D:\ECoG\timeseries_roi_213\';
Type = {'11';'22';'33';'12';'21';'23';'31'};

%{
for val_loop=156:200

Signal = [];State=[];
Nature = [];
dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(val_loop,'%04d')]);


for tl =1:numel(Type)
    x = load(fullfile(dest,['WO_Signal',Type{tl},'.mat']));
    Signal = cat(3,Signal,x.Y);
    x = load(fullfile(dest,['WO_State',Type{tl},'.mat']));
    State = cat(3,State,x.S);
    N = tl*ones(size(x.S,3),1);
    Nature = cat(1,Nature,N);
    if tl==6 | tl==7
        1;
    end
end
x=[];

X_RAW = {};
Y_RAW = {};
N = size(Signal,3);
for sl=1:N
    X_RAW{sl,1} = Signal(:,:,sl)';
    Y_RAW{sl,1} = Nature(sl)';
end
Signal = [];
YTest = categorical(Nature);

for TrainLength = 5

    val_loop
    
    filename = fullfile(dest,...
        ['Infomation0_GapLength_',num2str(00,'%02d'),...
        '_TrainLength_',num2str(TrainLength,'%02d'),...
        '_Layers_',num2str(1,'%02d'),...
        '_numHiddenUnits_',num2str(25,'%03d'),'.mat']);
    
    x=load(filename);
    Infomation = x.A; Y=[];x=[];
    
    SensitivityACC = zeros(numel(Infomation),7,100+1);
    for il=1:numel(Infomation)
        net = Infomation(il).net;
        miniBatchSize = 16*8;
        
        YPred_RAW = classify(net,X_RAW, ...
            'MiniBatchSize',miniBatchSize, ...
            'SequenceLength','longest');
        for natureType = 1:7
            bin=Nature==natureType;
            SensitivityACC(il,natureType,1) = sum(YPred_RAW(bin) == YTest(bin))./numel(YTest(bin));
        end
        for senloop= 1:100
            Xs_RAW = X_RAW;
            for nl=1:numel(Xs_RAW);Xs_RAW{nl}(senloop,:)=0;end
            YPred_s = classify(net,Xs_RAW, ...
                'MiniBatchSize',miniBatchSize, ...
                'SequenceLength','longest');
            for natureType = 1:7
                bin=Nature==natureType;
                SensitivityACC(il,natureType,senloop+1) = sum(YPred_s(bin) == YTest(bin))./numel(YTest(bin));
            end
            Xs_RAW = [];
        end
    end
    
    save(fullfile(dest,['SensitivityACC_TrainLength_',num2str(TrainLength),'.mat']),'SensitivityACC');
    
end
clear *RAW Y*
end
%}

%% PCA sensitivity
TrainLength=5;
for val_loop=1:500
    dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(val_loop,'%04d')]);
    load(fullfile(dest,['SensitivityACC_TrainLength_',num2str(TrainLength),'.mat']));
    if val_loop==1;Sacc=zeros(500,7,101);end
    Sacc(val_loop,:,:)=SensitivityACC;
end
Sacc =- Sacc+Sacc(:,:,1);
SensitivityACC = mean(Sacc,1);
yx=squeeze(SensitivityACC);
yx(:,1)=[];

ihdr = spm_vol(fullfile('H:\ECoG\pca','mask_RAM.nii'));
lmask = spm_read_vols(ihdr);
template = fullfile('H:\ECoG\Activation\','2nd','mean_T2.nii');
ihdr = spm_vol(template);
MaskFile = fullfile('H:\ECoG\Activation\','2nd','lmask_T2.nii');
cd('H:\ECoG\pca');
load('GroupWeight.mat')
filepath = fullfile(WholePath,'PCA_TEST','Sensitivity');
mkdir(filepath);
for ix=1:7
    Gcoeff = Gcoeff./std(Gcoeff,0,1);
    cd('H:\ECoG\pca');
    Sp = Gcoeff*yx(ix,:)'*100;
    fMRI = funmask(Sp,lmask);
    % sensitivity
    fMRI = fMRI/2+flip(fMRI,1)/2;
    hdr=ihdr; hdr.fname=fullfile(filepath,['Sensitivity_',num2str(ix),'.nii']);
    hdr.dt=[16 0];
    spm_write_vol(hdr,fMRI);
    
    % FDR mask
    y = squeeze(Sacc(:,ix,2:end));
    Spy = Gcoeff*y'*100;
    Vx = funmask(Spy,lmask);Vx = Vx/2+flip(Vx,1)/2; 
    Spy = fmask(Vx,lmask);
    [~,py] = ttest(Spy',0,'tail','right');
    [p_fdr, p_masked] = fdr( py, 10^-4,'nonParametric');
    d=zeros(1,size(Spy,1));
    for is=1:size(Spy,1); d(is) = computeCohen_d(Spy(is,:)',0);end
    Vmask = funmask((py<p_fdr & d>0.3)',lmask);
    %hdr=ihdr; hdr.fname=fullfile(filepath,['Mask_Sensitivity_',num2str(ix),'.nii']);
    %hdr.dt=[16 0];
    %spm_write_vol(hdr,Vmask);
    Cs = min(fMRI(Vmask==1));
    if isempty(Cs);Cs=5;end
    Cs
    
    
    barlim=[-5 -Cs Cs 5];
    %if ix==5|ix==4;barlim=[-1 -.01 0.01 1]*40;end
    %if ix==6|ix==7;barlim=[-1 -.01 0.01 1]*80;end
    
    Colormap(   'statfile',fullfile(filepath,['Sensitivity_',num2str(ix),'.nii']),...
        'bgmfile',template,...
        'slice',9:1:32,...
        'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
        'bar_value',barlim,...
        'dest',filepath,...
        'mapname',['Sensitivity_',num2str(ix)],...
        'cluster',20);

    barlim=[-5 -6 1.5 5];
    Colormap(   'statfile',fullfile(filepath,['Sensitivity_',num2str(ix),'.nii']),...
        'bgmfile',template,...
        'slice',9:2:34,...
        'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
        'bar_value',barlim,...
        'dest',filepath,...
        'mapname',['xSensitivity_',num2str(ix)],...
        'cluster',20);
    
    barlim=Cs;
    barmax = 5;
    Func3D=fullfile(filepath,['Sensitivity_',num2str(ix),'.nii']);
    I0 = spm_read_vols(spm_vol(Func3D));
    
    I = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Template_Mouse_X20.nii')));
    I0x = imresize3(I0,size(I),'linear');
    I0x = flip(flip(I0x,2),3);
    I0x(I0x<barlim)=barlim;
    I0x(I0x>barmax)=barmax;
    I0x(1:10,1:10,1:10)=barmax;
    
    for vloop=1:2
        config = struct('CameraPosition',[-4 0 0],...
                        'CameraUpVector',[0,0,1],...
                        'CameraTarget',[0 0 0],...
                        'CameraViewAngle',15,...
                        'BackgroundColor',[0 0 0]+1,...
                        'Renderer','MaximumIntensityProjection',... %%VolumeRendering
                        'Alphamap',(0:255)'/255,...
                        'Lighting',0,...
                        'IsosurfaceColor',[1 1 1]*0.5,...
                        'Isovalue',0.5);
        config.Colormap = [gray(256-20);ones(20,3)];
        if vloop==1
            config.CameraPosition=[0 -3 0];
            config.CameraUpVector = [-1 0 0];
        end
        if vloop==2
            config.CameraPosition=[-3 0 0];
            config.CameraUpVector = [0 1 0];
        end

        F1=figure;volshow(I,config);
        f1_rgb = getframe(F1);


        config.Colormap = MY_viridis_colormap(256)/256;
        config.Lighting = 0;
        F2=figure;volshow(I0x,config);
        f2_rgb = getframe(F2);

        f1 = f1_rgb.cdata;
        f2 = f2_rgb.cdata;
        f = f1*(1/2)+f2*(1/2);
        if vloop==1;F01=f;end
        if vloop==2;F02=f;end
    end
    close all
    Fa = uint8(255*ones(610,420+420,3));
    Fa(51:end,1:420,:) = flip(permute(F02,[2 1 3]),1);
    Fa(1:315,1+420:end,:) = imresize(F01,0.75);
    
    imwrite(uint8(Fa),fullfile(filepath,['Sensitivity_',num2str(ix),'_3D.tiff']));
    close all;
    
end



