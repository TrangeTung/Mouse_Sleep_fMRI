
clc;clear;
WholePath = 'F:\ECoG\m10_2\';
cd(WholePath)
hdr = spm_vol('2dseq.nii,1');
img = spm_read_vols(hdr);

lmask = spm_read_vols(spm_vol('EPI_mask.nii'));
for sl=22:-1:1
    %X = imread(['2dseq',num2str(sl,'%02d'),'.tif']);
   X = flip(img(:,:,sl)',1);
   Q = flip(lmask(:,:,sl)',1);
%    SE = strel('disk',3);
%    Bil = imdilate(Bil,SE);
%    Q = X>35; Q = Q.*(~Bil);
%    Q = imdilate(Q,SE);
%    SE = strel('disk',6);
%    Q = imerode(Q,SE);
%    SE = strel('disk',3);
%    Q = imdilate(Q,SE);
   
   y = repmat(Q/max(Q(:))*255,[1,1,3])/3;
   y(:,:,[2,3])=0;
   
   Y = uint8(X/3000*255+y);
   Y(:,[1:10,81:90],:)=[];
   imwrite(Y,['Cut_',num2str(sl),'.tiff']);
end



