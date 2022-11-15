clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
cd(codepath)
WholePath = 'H:\ECoG\ripple-average\';


Type = {'Awake SWR';'NREM SWR';...
    'spindle-uncoupled SWR';'spindle-coupled SWR';...
    'SWR-uncoupled spindle';'NREM spindle'};

v = VideoWriter(fullfile(WholePath,['Event_X1.avi']));
v.FrameRate=1;
open(v);


for f=1:13

    I = [];

    for tl = 1:numel(Type)

        if tl==1;    Func3D = fullfile(WholePath,'awake_ripple','awake_ripple.nii');end
        if tl==2;    Func3D = fullfile(WholePath,'nrem_ripple','nrem_ripple.nii');end
        if tl==4;    Func3D = fullfile(WholePath,'spindle_coupled','couple.nii');end
        if tl==3;    Func3D = fullfile(WholePath,'spindle_uncoupled','uncouple.nii');end
        if tl==5;    Func3D = fullfile(WholePath,'spindle','spindle.nii');end
        if tl==6;    Func3D = fullfile(WholePath,'nrem_spindle','nrem_spindle.nii');end


        close all
        Temp3D = fullfile(WholePath,'Template_Mouse_v38.nii');

        Img_RGB=Colormap(   'statfile',[Func3D,',',num2str(f)],...
                        'bgmfile',Temp3D,...
                        'slice',9:2:27,...
                        'bar_value',[-.1 -.02 .02 .1]*5,...
                        'denoi_profile',fullfile('H:\ECoG\pca','mask_RAM.nii'),...
                        'dest',fullfile(WholePath,'spindle'),...
                        'mapname',['s'],...
                        'cluster',10);
        Ix = uint8(Img_RGB);
        Iy = uint8(zeros([130 1140 3]));
        %if tl==4;
        Iy = insertText(Iy,[size(Ix,2)*0.4,0],['Time: ',num2str((f-7)*2),' s'],...
            'fontsize',25,'BoxColor','black','TextColor','yellow');
        %end
        Iy = insertText(Iy,[size(Ix,2)*0,0],Type{tl},...
            'fontsize',30,'BoxColor','black','TextColor','white');

        Iy(end-size(Ix,1)+1:end,:,:) = Ix;
        
        I = cat(1,I,Iy);
        
    end
    
    writeVideo(v,I);

    %{
    for f = 1:13
        
        FrameI = cell(1,1);
        
        
        
        Temp3D = fullfile(WholePath,'2nd','Template_Mouse_v38.nii');
        barlim = [-.101 .101]*5;
        %if sl==3|sl==4;barlim = [-.301 .301];end
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
        
        F0 = figure('color',[1 1 1]); axis off
        colormap('jet');caxis(barlim);
        c=colorbar('South','FontSize',18,'linewidth',1.5);
        c.Label.String = 'BOLD signal (z)';
        Ip = getframe(F0);Ip=Ip.cdata(201:end,:,:);
        I3(600+(1:size(Ip,1)),:,:)=Ip;

        I = cat(2,I0,I3);
        
        J = insertText(I,[0,size(I,2)*.0],['Time: ',num2str((f-7)*2),' s'],...
            'fontsize',30,'BoxColor','white');
        
        FrameI{1,1}=J;
        
        1;
        
        close all
        
        W = cell2mat(FrameI);
        Wr = W;%imresize(W,[1080 1920],'bilinear');

        writeVideo(v,Wr);

        
        
    end
    %}
    
    
    
end
    close(v);
