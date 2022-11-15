function [Fintra_dorsal,Fintra_lateral,Finter_dorsal,Finter_lateral] = ...
    MY_3D_node_CCmap(CCmatrix,Template_NII,Labels_NII,Labels_Excel,colorbar)


[~,~,CellData] = xlsread(Labels_Excel);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
T = spm_read_vols(spm_vol(Template_NII));
L = spm_read_vols(spm_vol(Labels_NII));



LOC=[];
for ls = 1:size(ExpTable,1)

    for sid = 1:2
        
        Mask0 = double(L==(sid-1)*size(ExpTable,1)+ls);
        Mask0 = smooth3(Mask0,'box',9);
        Mask0 = Mask0>0.2;

        Mask = Mask0;
        %if sid==1;Mask(1:114,:,:)=0;end
        %if sid==2;Mask(115:228,:,:)=0;end
        ind=find(Mask(:)==1);
        [I,J,K]=ind2sub(size(Mask),ind);
        Center = round([median(sort(I)),median(sort(J)),median(sort(K))]);
        LOC((sid-1)*size(ExpTable,1)+ls,:)=Center;
    end
end


%% intra-hemisphere
X = CCmatrix;
X(1:size(ExpTable,1),1+size(ExpTable,1):end)=nan;
X(1+size(ExpTable,1):end,1:size(ExpTable,1))=nan;

F = figure('color',[1 1 1]); %
ALL = T>191919;
p0 = patch(isosurface(ALL,0.5));
isonormals(ALL,p0);
p0.FaceColor = 'black';  p0.EdgeColor = 'none';
p0.FaceAlpha=0.02;     axis equal; axis tight
% axis off; 
hold on;

gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
ColorVec = [gdmap;flipud(drmap);];
lims = colorbar;
LineW = [128:-1:1,1:128]/16;

y = tril(X,-1);
for i=1:size(y,1)
    for j=i+1:size(y,2)
        value = y(j,i);
        if ~isnan(value) & value~=0
            
            value = y(j,i);
            pin = round ((value-lims(1)) / ((lims(2)-lims(1))/256));
            if pin>256;pin=256;end
            if pin<1;pin=1;end
            lis=LineW(pin);
            
            x_ = LOC(i,:);y_=LOC(j,:);
            h=plot3([x_(2),y_(2)],[x_(1),y_(1)],[x_(3),y_(3)],...
                'Color',ColorVec(pin,:),'linewidth',lis);
            h.Color(4)=0.3;
            hold on;
        end
    end
end
axis off;


set(F,'Position', [680 87 728 891])
view([-90 0]); f_ = getframe(F);
Fintra_dorsal = f_.cdata;
set(F,'Position', [680 87 525 891])
view([-180 0]);f_ = getframe(F);
Fintra_lateral = f_.cdata;



%% inter-hemisphere
X = CCmatrix;
X(1:size(ExpTable,1),1:size(ExpTable,1))=nan;
X(1+size(ExpTable,1):end,1+size(ExpTable,1):end)=nan;

F = figure('color',[1 1 1]); %
ALL = T>191919;
p0 = patch(isosurface(ALL,0.5));
isonormals(ALL,p0);
p0.FaceColor = 'black';  p0.EdgeColor = 'none';
p0.FaceAlpha=0.02;     axis equal; axis tight
% axis off; 
hold on;

gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
ColorVec = [gdmap;flipud(drmap);];
lims = colorbar;
LineW = [128:-1:1,1:128]/16;

y = tril(X,-1);
for i=1:size(y,1)
    for j=i+1:size(y,2)
        value = y(j,i);
        if ~isnan(value) & value~=0
            
            value = y(j,i);
            pin = round ((value-lims(1)) / ((lims(2)-lims(1))/256));
            if pin>256;pin=256;end
            if pin<1;pin=1;end
            lis=LineW(pin);
            
            x_ = LOC(i,:);y_=LOC(j,:);
            h=plot3([x_(2),y_(2)],[x_(1),y_(1)],[x_(3),y_(3)],...
                'Color',ColorVec(pin,:),'linewidth',lis);
            h.Color(4)=0.3;
            hold on;
        end
    end
end
axis off;


set(F,'Position', [680 87 728 891])
view([-90 0]); f_ = getframe(F);
Finter_dorsal = f_.cdata;
set(F,'Position', [680 87 525 891])
view([-180 0]);f_ = getframe(F);
Finter_lateral = f_.cdata;


end