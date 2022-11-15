clc;clear
WholePath = 'F:\ECoG\pc_yyl';


type={'';'_rp_6';'_rp_12';'_1';'_10';'_20';'_30';'_40'};
colorVec = {[0,0,0];[164,12,94];[146,8,131];[203,224,214];[195,195,174];[102,166,134];[42,130,86];[0,101,48]};
%% DVARS pdf
%Ds = zeros(7200*46,numel(type));
for idx=1:46
    cd(fullfile(WholePath,num2str(idx,'%02d')));
    
    rp = load('rp_sm2dseq.txt');
    drp = rp - [rp(1,:);rp(1:end-1,:)];
    drp(1:3)=drp(1:3)/20;
    drp(4:6)=drp(4:6)*5;
    FD=sqrt(drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2);
    
    for tl=1:numel(type)
        x = load(['G_S',type{tl},'.txt']);
        if tl==1;x=rp(:,3);end
        [pxx,f,~] = periodogram(x,rectwin(length(x)),length(x),1/2);
                
        psd(idx,tl,:) = pxx;
        
    end
    
    
  
end

F = figure;
for tl=1:numel(type)
    pxx = squeeze(mean(psd(:,tl,:),1));
    semilogx(f,10*log10(pxx))
    hold on;
end
xlim([0 max(f)])
xlabel('Hz')
ylabel('dB/Hz')
set(gca,'box','off','linewidth',.5,'tickdir','out');





clc;clear
WholePath = 'F:\ECoG\shang_m2_3';
filename = fullfile(WholePath,'snrsm2dseq.nii');
hdr = spm_vol(filename);
fMRI_4D = spm_read_vols_4D(hdr);

lmask = spm_read_vols_4D(spm_vol(fullfile(WholePath,'mask.nii')));

rp = load(fullfile(WholePath,'rp_sm2dseq.txt'));
drp = rp - [rp(1,:);rp(1:end-1,:)];
drp(1:3)=drp(1:3)/20;
drp(4:6)=drp(4:6)*5;
FD=sqrt(drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2);
PCs = load(fullfile(WholePath,'PC40.txt'));

RegBasFunc1 = [ones([length(rp),1]) rp];
RegBasFunc2 = [ones([length(rp),1]) rp drp ];
RegBasFunc3 = [ones([length(rp),1]) rp drp PCs(:,1)];
RegBasFunc4 = [ones([length(rp),1]) rp drp PCs(:,1:10)];
RegBasFunc5 = [ones([length(rp),1]) rp drp PCs(:,1:20)];
RegBasFunc6 = [ones([length(rp),1]) rp drp PCs(:,1:30)];
RegBasFunc7 = [ones([length(rp),1]) rp drp PCs(:,1:40)];


data_ready_regress0 = fmask(fMRI_4D,lmask);   
data_ready_regress1 = data_ready_regress0*0;
data_ready_regress2 = data_ready_regress0*0;
data_ready_regress3 = data_ready_regress0*0;
data_ready_regress4 = data_ready_regress0*0;
data_ready_regress5 = data_ready_regress0*0;
data_ready_regress6 = data_ready_regress0*0;
data_ready_regress7 = data_ready_regress0*0;

for iii = 1:size(data_ready_regress0,1)
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc1);
    data_ready_regress1(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc2);
    data_ready_regress2(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc3);
    data_ready_regress3(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc4);
    data_ready_regress4(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc5);
    data_ready_regress5(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc6);
    data_ready_regress6(iii,:) = Residual + Beta(1);
    [Beta,~,Residual] = regress(double(data_ready_regress0(iii,:))',RegBasFunc7);
    data_ready_regress7(iii,:) = Residual + Beta(1);
end

[CC0,p0] = corr(data_ready_regress0',FD(:));
[CC1,p1] = corr(data_ready_regress1',FD(:));
[CC2,p2] = corr(data_ready_regress2',FD(:));
[CC3,p3] = corr(data_ready_regress3',FD(:));
[CC4,p4] = corr(data_ready_regress4',FD(:));
[CC5,p5] = corr(data_ready_regress5',FD(:));
[CC6,p6] = corr(data_ready_regress6',FD(:));
[CC7,p7] = corr(data_ready_regress7',FD(:));


cd(WholePath); clear data_*

C3D = funmask(CC0,lmask);p3D = funmask(p0,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC0.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p0.nii';spm_write_vol(ihdr,p3D);

C3D = funmask(CC1,lmask);p3D = funmask(p1,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC1.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p1.nii';spm_write_vol(ihdr,p3D);

C3D = funmask(CC2,lmask);p3D = funmask(p2,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC2.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p2.nii';spm_write_vol(ihdr,p3D);

C3D = funmask(CC3,lmask);p3D = funmask(p3,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC3.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p3.nii';spm_write_vol(ihdr,p3D);

C3D = funmask(CC4,lmask);p3D = funmask(p4,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC4.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p4.nii';spm_write_vol(ihdr,p3D);

C3D = funmask(CC5,lmask);p3D = funmask(p5,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC5.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p5.nii';spm_write_vol(ihdr,p3D);

C3D = funmask(CC6,lmask);p3D = funmask(p6,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC6.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p6.nii';spm_write_vol(ihdr,p3D);

C3D = funmask(CC7,lmask);p3D = funmask(p7,lmask);
ihdr=hdr(1);ihdr.dt=[16,0];
ihdr.fname='CC7.nii';spm_write_vol(ihdr,C3D);
ihdr.fname='p7.nii';spm_write_vol(ihdr,p3D);

type={'';'_rp_6';'_rp_12';'_1';'_10';'_20';'_30';'_40'};
for loop=0:7
    cd(WholePath);
    img = spm_read_vols(spm_vol(['CC',num2str(loop),'.nii']));
    pmask = spm_read_vols(spm_vol(['p',num2str(loop),'.nii']));
    
    slice = [11:5:30];
    bar_value = [0.1 0.3];
    defaultMap = gray(64);%
    
    Func_Img_3D = abs(img);
    map_nothre_reshape = reshape(permute(Func_Img_3D(:,:,slice),[2 1 3]),[size(Func_Img_3D,2) size(Func_Img_3D,1)*numel(slice)]);

    nothrebar = [min(bar_value(:)) max(bar_value(:))];
    tmap = ones([size(map_nothre_reshape),3]);
    map_normalize = (map_nothre_reshape-nothrebar(1))/(nothrebar(2)-nothrebar(1));
    map_normalize(map_normalize<=0) = 0.000000001;
    map_normalize(map_normalize>=1) = 1;
    map_normalize(isnan(map_normalize)) = 0.000000001;
    tmap(:,:,1) = reshape(defaultMap(ceil(map_normalize*64),1)*255,size(map_normalize));
    tmap(:,:,2) = reshape(defaultMap(ceil(map_normalize*64),2)*255,size(map_normalize));
    tmap(:,:,3) = reshape(defaultMap(ceil(map_normalize*64),3)*255,size(map_normalize));

    mask_reslice = ~lmask(:,:,slice);
    mask_reshape = reshape(permute(mask_reslice,[2 1 3]),[size(mask_reslice,2) size(mask_reslice,1)*numel(slice)]);
    mask_pixel = (mask_reshape == 1);
    mask_index = find(double(mask_pixel)==1);
    tmap(mask_index+numel(mask_reshape)*0) = 0;
    tmap(mask_index+numel(mask_reshape)*1) = 0;
    tmap(mask_index+numel(mask_reshape)*2) = 0;
    tmap = flip(tmap,1);
    
    for sl=1:numel(slice)
        col=sl*size(tmap,2)/numel(slice)+(-size(tmap,2)/numel(slice):-1)+1;
        if sl==1
            I=tmap(:,col,:);
        else
            I=cat(1,I,tmap(:,col,:));
        end
    end
    
    filename = ['Map_',type{loop+1},'.tiff'];
    imwrite(uint8(I),filename);
    
end

