clc;clear
WholePath = 'F:\ECoG\pc_yyl';


type={'';'_rp_6';'_rp_12';'_1';'_10';'_20';'_30';'_40'};
colorVec = {[0,0,0];[164,12,94];[146,8,131];[203,224,214];[195,195,174];[102,166,134];[42,130,86];[0,101,48]};
%% DVARS pdf
%Ds = zeros(7200*46,numel(type));
CC = zeros(46,numel(type));
for idx=1:46
    cd(fullfile(WholePath,num2str(idx,'%02d')));
    
    rp = load('rp_sm2dseq.txt');
    drp = rp - [rp(1,:);rp(1:end-1,:)];
    drp(1:3)=drp(1:3)/20;
    drp(4:6)=drp(4:6)*5;
    FD=sqrt(drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2);
    
    for tl=1:numel(type)
        dvars = load(['dvars',type{tl},'.txt']);
        CC(idx,tl) = corr(dvars(:),FD(:));
        %Ds((idx-1)*numel(dvars)+(1:numel(dvars)),tl) = dvars;
        
          shift = [-5:5];
            for sl=1:11
                X = circshift(dvars(:),shift(sl),1);
                CCs(idx,tl,sl) = corr(FD,X);
            end
        
    end
    
    
  
end

F = figure;
for tl=1:numel(type)
    histogram(CC(:,tl),20,'Normalization','count','LineWidth',.5,...
        'DisplayStyle','stairs','edgecolor',colorVec{tl}/255);
    hold on;
end
xlim([0 1])
set(gca,'box','off','linewidth',.5,'tickdir','out');
