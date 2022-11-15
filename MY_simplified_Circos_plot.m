function F=MY_simplified_Circos_plot(CCMatrix,AtlasTable,colorbar)

opt.label = [cell(size(flip(AtlasTable.abbre,1)));AtlasTable.abbre;];
gdmap = [(0:9)'/9,(0:9)'/9,ones(10,1)];
drmap = [ones(10,1),(0:9)'/9,(0:9)'/9];
defaultMap = [gdmap;flipud(drmap);];
opt.ccolor = defaultMap;
CV = cat(1,AtlasTable.ColorVec);
CV_ = cellfun(@(x) str2num(x), CV, 'UniformOutput', false);
opt.ncolor = [flipud(cell2mat(CV_));cell2mat(CV_)]/255;
MajorRegion = cat(1,AtlasTable.Major_Region);

Labs = {'Isocortex_primary';'Isocortex_HigherOrder';'HPF';'STR';'PAL';'TH';'HY';'MB'};
Defs = {[1:18];[18:31];[32:36];[38:41];[42:45];[46:55];[56:59];[60:69]};
opt.label = Labs;

%% 
Y = CCMatrix; 

MatIntra = nan([numel(Labs),numel(Labs)]);
for sl=1:numel(Labs)
    for tl=1:numel(Labs)
        mC = nanmean(nanmean(Y(Defs{sl},Defs{tl})));
        MatIntra(sl,tl) = mC;
    end
end

MatInter = nan([numel(Labs),numel(Labs)]);
for sl=1:numel(Labs)
    for tl=1:numel(Labs)
        mC = nanmean(nanmean(Y(Defs{sl},Defs{tl}+numel(MajorRegion))));
        MatInter(sl,tl) = mC;
    end
end
XZ = (((MatInter)));
YZ = (fliplr(tril(MatIntra)));

Mat = cat(1,XZ);
DirecMat = cat(2,Mat);
Y = DirecMat;

CV = {[39,143,87];[39,143,87]*.7;[123,75,30];[238,173,77];[183,28,37];[220,0,0];[153,0,255];[154,27,91]};
opt.ncolor = [cell2mat(CV)]/255;


Y_=Y;Y_(Y_<0)=0; Y_pos = Y_/abs(max(colorbar)); Y_pos(Y_pos>1)=1; Y_pos(isnan(Y_pos))=0;
Y_=Y;Y_(Y_>0)=0; Y_neg = Y_/abs(min(colorbar)); Y_neg(Y_neg<-1)=-1;Y_neg(isnan(Y_neg))=0;
Yf = Y_pos+Y_neg;
F = MY_Circos_plot(Yf,opt);

% close all

end