clc;clear
close all
codepath = 'F:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'F:\ECoG\pca\';
cd(codepath);

ALL=[];
Seg_aw=[];Seg_nrem=[];Seg_rem=[];
for idx=1:46
    
    
    load(fullfile('F:\ECoG\state_check\',['s',num2str(idx),'.mat']));
    ALL = all_state;
    
    x = ALL==1; y = x - [x(1);x(1:end-1)]; 
    L1 = find(y==1);
    L2 = find(y==-1);
    if L1(1)>L2(1);L1=[1;L1];end
    NUM = min([numel(L1),numel(L2)]);
    Seg = L2(1:NUM)-L1(1:NUM);
    
    if find(Seg>3000)
        figure;plot(ALL)
       1; 
    end
    
    Seg_aw = cat(1,Seg_aw,Seg);

    x = ALL==3; y = x - [x(1);x(1:end-1)]; 
    L1 = find(y==1);
    L2 = find(y==-1);
    if L1(1)>L2(1);L1=[1;L1];end
    NUM = min([numel(L1),numel(L2)]);
    Seg = L2(1:NUM)-L1(1:NUM);
    Seg_nrem = cat(1,Seg_nrem,Seg);

    x = ALL==4; y = x - [x(1);x(1:end-1)]; 
    L1 = find(y==1);
    L2 = find(y==-1);
    if ~isempty(L1)
        if L1(1)>L2(1);L1=[1;L1];end
        NUM = min([numel(L1),numel(L2)]);
        Seg = L2(1:NUM)-L1(1:NUM);
        Seg_rem = cat(1,Seg_rem,Seg);
    end
    
    
    
end


figure;
histogram(Seg_aw,'BinWidth',5,'DisplayStyle','stairs');hold on;
histogram(Seg_nrem,'BinWidth',5,'DisplayStyle','stairs');hold on;
histogram(Seg_rem,'BinWidth',5,'DisplayStyle','stairs');hold on;
xlim([0 500])



figure;
plot(sort(Seg_aw));hold on;
plot(sort(Seg_nrem));hold on;
plot(sort(Seg_rem));hold on;
set(gca,'yScale','log');

