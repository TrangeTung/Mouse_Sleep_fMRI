%% temporal denpendency
%
clc;clear
addpath(genpath('H:\ECoG\code_Trange'));
WholePath = 'H:\ECoG\timeseries_roi_213\';
Type = {'11';'22';'33';'12';'21';'23';'31'};

LayerNumber = 1;
numHiddenUnits = 25;
ACCall = [];
ACCall_shuffle = [];
%{
for shf_loop=1:500
    
shf_loop
    %% actual
    Signal = [];State=[];
    Nature = [];
    dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(shf_loop,'%04d')]);

    for tl =1:numel(Type)
        x = load(fullfile(dest,['WO_Signal',Type{tl},'.mat']));
        Signal = cat(3,Signal,x.Y);
        x = load(fullfile(dest,['WO_State',Type{tl},'.mat']));
        State = cat(3,State,x.S);
        N = tl*ones(size(x.S,3),1);
        Nature = cat(1,Nature,N);
        if tl==6 | tl==7
            1;
        end
    end
    x=[];

    X_RAW = {};
    Y_RAW = {};
    N = size(Signal,3);
    for sl=1:N
        X_RAW{sl,1} = Signal(:,:,sl)';
        Y_RAW{sl,1} = Nature(sl)';
    end
    Signal = [];
    YTest = categorical(Nature);
    
    
for TrainLength = 5
for GapLength = 00:30
    cd(dest);
    
    filename = fullfile(dest,...
        ['Infomation0_GapLength_',num2str(00,'%02d'),...
        '_TrainLength_',num2str(TrainLength,'%02d'),...
        '_Layers_',num2str(LayerNumber,'%02d'),...
        '_numHiddenUnits_',num2str(numHiddenUnits,'%03d'),'.mat']);
    load(filename);
    
    Infomation = A; Y=[];
    
    %
    for il=1:numel(Infomation)
        net = Infomation(il).net;
        miniBatchSize = 64;
        
        RdN = Infomation(il).randomY;
        predicted_Time = 62;
        input_Time = 61+(-TrainLength+1:0)-GapLength;
        for sl=1:numel(RdN)
            XTrain{sl,1} = X_RAW{RdN(sl),1}(:,input_Time);
            YTrain(sl,1) = Y_RAW{RdN(sl),1}(1);
        end
        YTest = categorical(YTrain);
        
        YPred = classify(net,XTrain, ...
                'MiniBatchSize',miniBatchSize, ...
                'SequenceLength','longest');
        bin=Nature(RdN)==1;acc(il,1) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==2;acc(il,2) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==3;acc(il,3) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==4;acc(il,4) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==5;acc(il,5) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==6;acc(il,6) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==7;acc(il,7) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        
    end
    ACCall(TrainLength,GapLength+1,shf_loop,:) = acc;%mean(acc,1);
    %
end
end

%% shuffle
    Signal = [];State=[];
    Nature = [];
    dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(shf_loop,'%04d')]);

    for tl =1:numel(Type)
        x = load(fullfile(dest,['WO_Signal',Type{tl},'.mat']));
        Signal = cat(3,Signal,x.Y);
        x = load(fullfile(dest,['WO_State',Type{tl},'.mat']));
        State = cat(3,State,x.S);
        N = tl*ones(size(x.S,3),1);
        Nature = cat(1,Nature,N);
        if tl==6 | tl==7
            1;
        end
    end
    Nature = Nature(randperm(numel(Nature)));
    x=[];

    X_RAW = {};
    Y_RAW = {};
    N = size(Signal,3);
    for sl=1:N
        X_RAW{sl,1} = Signal(:,:,sl)';
        Y_RAW{sl,1} = Nature(sl)';
    end
    Signal = [];
    YTest = categorical(Nature);
    
    
for TrainLength = 5
for GapLength = 00:30
    cd(dest);
    
    filename = fullfile(dest,...
        ['Infomation0_GapLength_',num2str(00,'%02d'),...
        '_TrainLength_',num2str(TrainLength,'%02d'),...
        '_Layers_',num2str(LayerNumber,'%02d'),...
        '_numHiddenUnits_',num2str(numHiddenUnits,'%03d'),'.mat']);
    load(filename);
    
    Infomation = A; Y=[];
    
    %
    for il=1:numel(Infomation)
        net = Infomation(il).net;
        miniBatchSize = 64;
        
        RdN = Infomation(il).randomY;
        predicted_Time = 62;
        input_Time = 61+(-TrainLength+1:0)-GapLength;
        for sl=1:numel(RdN)
            XTrain{sl,1} = X_RAW{RdN(sl),1}(:,input_Time);
            YTrain(sl,1) = Y_RAW{RdN(sl),1}(1);
        end
        YTest = categorical(YTrain);
        
        YPred = classify(net,XTrain, ...
                'MiniBatchSize',miniBatchSize, ...
                'SequenceLength','longest');
        bin=Nature(RdN)==1;acc(il,1) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==2;acc(il,2) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==3;acc(il,3) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==4;acc(il,4) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==5;acc(il,5) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==6;acc(il,6) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        bin=Nature(RdN)==7;acc(il,7) = sum(YPred(bin) == YTest(bin))./numel(YTest(bin));
        
    end
    ACCall_shuffle(TrainLength,GapLength+1,shf_loop,:) = acc;%mean(acc,1);
    %
end
end


end

cd(fullfile(WholePath,'PCA_TEST'));
save('ACC.mat','ACC*');
%}

cd(fullfile(WholePath,'PCA_TEST'));
load('ACC.mat');
close all
for bin = 4:7
    batch_length = 5;
    X = squeeze(ACCall(batch_length,:,:,bin)); 
    Y = squeeze(ACCall_shuffle(batch_length,:,:,bin)); 
% % % % % % % % % % % % % % % % % %     if bin==7; Y=Y*7/2;end
% % % % % % % % % % % % % % % % % %     if bin==5 | bin==6; Y=Y*7/3;end
% % % % % % % % % % % % % % % % % %     if bin==4; Y=Y*7/3.5;end
    
    F=figure('Position', [138*(bin-4)*3 221 382 757*3/4]);
    subplot(3,1,1:2);
    x0 = (1:size(X,1))*2-2;
    Ms = median(X,2); SEMs = std(X,0,2)*2;%/sqrt(size(X,2));
    plot(x0,Ms,'r-','linewidth',2); hold on;
%     patch('XData',[x0,fliplr(x0)],'YData',[Ms+SEMs;flipud(Ms-SEMs)],...
%         'EdgeColor','none','FaceColor','r','FaceAlpha',0.25);
    x0 = (1:size(Y,1))*2-2;
    Ms = mean(Y,2); SEMs = std(Y,0,2)*2;%/sqrt(size(Y,2));
    plot(x0,Ms,'K-','linewidth',2); hold on;
    patch('XData',[x0,fliplr(x0)],'YData',[Ms+SEMs;flipud(Ms-SEMs)],...
        'EdgeColor','none','FaceColor','K','FaceAlpha',0.25);
    set(gca,'XDir','reverse');
    xlabel('Gap (s)'); ylabel('Accuracy (%)');
    xlim([-2 30+2]);ylim([0 1])
    set(gca,'box','off','linewidth',1,'tickdir','out','fontsize',12);
    set(gca,'yaxislocation','right','ticklength',[0.04 1]);
    subplot(3,1,3);
    TPP = zeros(500,1);
    for tm=1:500
        T = max(find(X(:,tm)>Ms+SEMs));
        if isempty(T);T=nan;end
        TPP(tm) = (T-1)*2;
%        [~,p(tm)]=ttest(X(tm,:),Y(tm,:),'tail','right'); 
    end
    if bin==7;TPP(TPP>10)=[];end
    TPP(TPP>nanmedian(TPP)+nanstd(TPP)*2|TPP<nanmedian(TPP)-nanstd(TPP)*2)=nan;
    histogram(TPP,'BinWidth',2)
    %raincloud_plot(TPP,0,'box_on',1);
    set(gca,'XDir','reverse');
    xlabel('Gap (s)');
    xlim([-2 30+2]);
    set(gca,'box','off','linewidth',1,'tickdir','out','fontsize',12);
    set(gca,'yaxislocation','right');
    nanmean(TPP)
    saveas(F,fullfile(WholePath,'PCA_TEST',[Type{bin},'_acc.emf']));
end



1;
