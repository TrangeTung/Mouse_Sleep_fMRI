function [Fintra,Finter]=MY_modified_Circos_plot(CCMatrix,AtlasTable,colorbar)

opt.label = [cell(size(flip(AtlasTable.abbre,1)));AtlasTable.abbre;];
gdmap = [(0:9)'/9,(0:9)'/9,ones(10,1)];
drmap = [ones(10,1),(0:9)'/9,(0:9)'/9];
defaultMap = [gdmap;flipud(drmap);];
opt.ccolor = defaultMap;
CV = cat(1,AtlasTable.ColorVec);
CV_ = cellfun(@(x) str2num(x), CV, 'UniformOutput', false);
opt.ncolor = [flipud(cell2mat(CV_));cell2mat(CV_)]/255;
%% intra-hemisphere
Y = CCMatrix; 
Y(1:size(AtlasTable,1),size(AtlasTable,1)+1:end)=nan;
Y(size(AtlasTable,1)+1:end,1:size(AtlasTable,1))=nan;
Y(1:size(AtlasTable,1),1:size(AtlasTable,1)) = ...
    flip(flip(Y(1:size(AtlasTable,1),1:size(AtlasTable,1)),1),2);
Y_=Y;Y_(Y_<0)=0; Y_pos = Y_/abs(max(colorbar)); Y_pos(Y_pos>1)=1;
Y_=Y;Y_(Y_>0)=0; Y_neg = Y_/abs(min(colorbar)); Y_neg(Y_neg<-1)=-1;
Y_intra = nan(size(Y));
Y_intra(1:size(AtlasTable,1),1:size(AtlasTable,1))=Y_neg(1:size(AtlasTable,1),1:size(AtlasTable,1));
Y_intra(1+size(AtlasTable,1):end,1+size(AtlasTable,1):end)=flip(flip(Y_pos(1:size(AtlasTable,1),1:size(AtlasTable,1)),1),2);
Fintra = MY_Circos_plot(Y_intra,opt);
%% inter-hemisphere
Y = CCMatrix; 
Y(1:size(AtlasTable,1),1:size(AtlasTable,1))=nan;
Y(size(AtlasTable,1)+1:end,size(AtlasTable,1)+1:end)=nan;
Y(1:size(AtlasTable,1),size(AtlasTable,1)+1:end) = ...
    flip(Y(1:size(AtlasTable,1),size(AtlasTable,1)+1:end),1);
Y(size(AtlasTable,1)+1:end,1:size(AtlasTable,1)) = ...
    flip(Y(size(AtlasTable,1)+1:end,1:size(AtlasTable,1)),2);
Y_=Y;Y_(Y_<0)=0; Y_pos = Y_/abs(max(colorbar));Y_pos(Y_pos>1)=1;
Y_=Y;Y_(Y_>0)=0; Y_neg = Y_/abs(min(colorbar));Y_neg(Y_neg<-1)=-1;
Y_inter=Y_neg+Y_pos;
Finter = MY_Circos_plot(Y_inter,opt);

% close all

end