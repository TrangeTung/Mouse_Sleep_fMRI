clc;clear
codepath = 'F:\ECoG\';
path = 'F:\ECoG\Activation\';

X = spm_read_vols(spm_vol(fullfile(codepath,'Colormap_3Dviewer','Label_Mouse_X20.nii')));
Y = flip(permute(X,[3 2 1]),1);
Labels = Y;%imresize3(Y,size(Y)*2,'nearest');
fimgx = Labels*0;

for loop=1
    config = struct('CameraPosition',[-4 0 0],...
        'CameraUpVector',[0,0,1],...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',15,...
        'BackgroundColor',[0 0 0]+1);
    label = floor(Labels/6);
    if loop==1;Cs=[0 127];end
    if loop==2;Cs=[128 255];end
    if loop==3;Cs=[256 378];end
    
    label(label<Cs(1))=0;
    label(label>Cs(2))=0;
    Ns=Labels*0;Ns(:,:,size(Ns,3)/2+1:end)=1;
    label(Ns&Labels<=378)=0;
    label=round(label);
    label(label==92)=0;
    
    
    config.CameraPosition=[-3.45,-3.00,4.21]/2;
    config.CameraUpVector = [-1 0 0];
    F2=figure;labelvolshow(label,fimgx,config);
    f2_rgb = getframe(F2);
    imwrite(uint8(f2_rgb.cdata),fullfile(path,['Mask2X_',num2str(loop,'%02d'),'.tiff']));
    
    config.CameraPosition=[0 0 -3];
    config.CameraUpVector = [-1 0 0];
    F1=figure;labelvolshow(label,fimgx,config);
    f1_rgb = getframe(F1);
    imwrite(uint8(f1_rgb.cdata),fullfile(path,['Mask1X_',num2str(loop,'%02d'),'.tiff']));
    
end

label = Labels; label(Labels<640)=0;
label(Labels>872)=0;
label = label-640; label(label<0)=0;
label = round(label/3);
config.CameraPosition=[-3.45,-3.00,4.21]/2;
config.CameraUpVector = [-1 0 0];
F3=figure;labelvolshow(label,fimgx,config);
f3_rgb = getframe(F3);
imwrite(uint8(f3_rgb.cdata),fullfile(path,'Mask3.tiff'));

