function Fcdata = MY_Global_Connectivity_map(CCmatrix,Template_NII,Labels_NII,Labels_Excel,colorbar)

[~,~,CellData] = xlsread(Labels_Excel);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
MajorRegion = cat(1,ExpTable.Major_Region);
T = spm_read_vols(spm_vol(Template_NII));
L = spm_read_vols(spm_vol(Labels_NII));

X = CCmatrix;
GBC = X;
% GBC = ( nansum(X,2))/size(X,1);
% GBC = GBC/2+([GBC(numel(MajorRegion)+1:end);GBC(1:numel(MajorRegion))])/2;

gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
ColorVec = [gdmap;flipud(drmap);];
ColorVec = MY_viridis_colormap(256)/255;
lims = colorbar;

for gloop=1:3
    
    if gloop==1;Ss={'Isocortex'};end
    if gloop==2;Ss={'Hippocampal Formation';'Striatum';};end
    if gloop==3;Ss={'Thalamus';'Pallidum';'Hypothalamus';'Midbrain'};end
    
    %% figure
    F = figure('color',[1 1 1]); %
    subplot('position',[0.02,0.02,1-0.04,1-0.04]);
    ALL = T>191919;
    p = patch(isosurface(ALL,0.5));
    isonormals(ALL,p);
    p.FaceColor = 'black';  p.EdgeColor = 'none';
    p.FaceAlpha=0.02;     axis equal; axis tight
    % axis off;
    hold on;
    for ls = 1:numel(GBC)
        
        lss=ls;
        if ls>numel(MajorRegion);lss=ls-numel(MajorRegion);end
        switch MajorRegion{lss}
            case Ss
                
                value = GBC(ls);
                pin = round ((value-lims(1)) / ((lims(2)-lims(1))/256));
                if pin>256;pin=256;end
                if pin<1;pin=1;end
                
                Mask0 = L==ls;
                
                Mask0 = smooth3(Mask0,'box',11);
                L0 = Mask0>0.2;
                
                ALL = double(L0);
                p0 = patch(isosurface(ALL,0.5));
                isonormals(ALL,p0);
                p0.FaceColor = ColorVec(pin,:);
                p0.EdgeColor = 'none';
                p0.FaceAlpha=.25;
                axis equal; axis tight
                % axis off;
                hold on;
        end
        
    end
    set(F,'Position', [680 87 728 891])
    view([-90 0]);
    if gloop~=1;view([-180 0]);end
    axis off;
    
    f_ = getframe(F);
    Fcdata{:,gloop} = f_.cdata;
end

end