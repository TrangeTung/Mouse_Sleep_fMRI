%% 3D fitting x:gap length y:train length z:accuracy/loss
clear;clc;
% close all
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'D:\ECoG\timeseries_roi_213\';
Type = {'11';'22';'33';'12';'21';'23';'31'};
%{
for val_loop=1:100
Signal = [];State=[];
Nature = [];
dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(val_loop,'%04d')]);

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

TrainLength = 5;
GapLength = 0;
for LayerNumber = 1:3
    count=0;
    
    for numHiddenUnits = 25:25:200
        
        filename = fullfile(dest,...
            ['Infomation0_GapLength_',num2str(GapLength,'%02d'),...
            '_TrainLength_',num2str(TrainLength,'%02d'),...
            '_Layers_',num2str(LayerNumber,'%02d'),...
            '_numHiddenUnits_',num2str(numHiddenUnits,'%03d'),'.mat']);
        load(filename);
        Infomation = A;
        
        
        for il=1:numel(Infomation)
            
            %RdN = Infomation(il).randomY;
            RdN = randperm(numel(X_RAW),2000);
%             predicted_Time = 62;
%             input_Time = 61+(-TrainLength+1:0)-GapLength;
%             for sl=1:numel(RdN)
%                 XTrain{sl,1} = X_RAW{RdN(sl),1}(:,input_Time);
%                 YTrain(sl,1) = Y_RAW{RdN(sl),1}(1);
%             end
%             YTest = categorical(YTrain);
%             
%             
%             net = Infomation(1).net;
%             YPred = classify(net,XTrain, ...
%                 'MiniBatchSize',16*8, ...
%                 'SequenceLength','longest');
%             
%             acc = sum(YPred == YTest)./numel(YTest);
            acc=Infomation.TestACC;
            ACC(il) = acc;
        end
        %         for il=1:numel(Infomation);a(il,:) = Infomation(il).ValidationAccuracy;end
        %         ACC = max(a);
        %         for il=1:numel(Infomation);ACC(il) = Infomation(il).TestACC;end
        pin = min(find(ACC==max(ACC)));
        for il=1:numel(Infomation);Y(:,il) = Infomation(il).ValidationLoss;end
        Y = nanmean(Y,2);
        
        count=count+1;
        ACCmax(LayerNumber,count,val_loop) = mean(ACC);
        Lossmin(LayerNumber,count,val_loop) = min(Y(end));
        1;
    end
end
end
numHiddenUnits = [25:25:200];
F = figure('Position', [680 599 419 379]);
for loop=1:3
X = squeeze(ACCmax(loop,:,:));
Ms = mean(X,2);
SEMs = std(X,0,2)/sqrt(size(X,2));
errorbar(numHiddenUnits,Ms,SEMs,'linewidth',1);
hold on
end
ylabel('Accuracy');
ylim([0.9 0.98]);
xlim([0 225])
xlabel('No. of Hidden Units');grid on;
set(gca,'box','off','tickdir','out','linewidth',2,'fontsize',15);





%}


%{
for val_loop=1:100
Signal = [];State=[];
Nature = [];
dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(val_loop,'%04d')]);

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

LayerNumber = 1;
numHiddenUnits = 25;

for TrainLength = 5:15
for GapLength = 00:15
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
        miniBatchSize = 16*8;
        
        RdN = 1:numel(X_RAW);%Infomation(il).randomY;
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
        
        for pr=1:7
            for te=1:7
                N_ = numel(find(double(YPred)==pr & double(YTest)==te))...
                    /numel(find(double(YTest)==te));
                Ctest(pr,te,il) = N_;
            end
        end
        
    end
    ACCall(val_loop,TrainLength-4,GapLength+1,:,:) = acc;%mean(acc,1);
    %
    ConfMat(val_loop,TrainLength-4,GapLength+1,:,:) = mean(Ctest,3);
end
end
end

Z = squeeze(mean(ConfMat(:,1,1,:,:),1))*100;
B = squeeze(mean(mean(ACCall,5),1));
[X,Y] = meshgrid(2*(0:15),2*(5:15));
B(B==0)=nan;
zData = B(~isnan(B));
xData = X(~isnan(B));
yData = Y(~isnan(B));

ft = fittype( 'poly55' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

F = figure('Position', [680 621 376 357]);
h = plot( fitresult, [xData, yData], zData );
xlabel('Gap (s)');ylabel('Train length (s)');
zlabel('Train accuracy (%)');
h(1).EdgeColor = 'none';
h(1).FaceAlpha = 0.5;
h(1).FaceColor = 'interp';
h(2).MarkerSize = 6;
h(2).MarkerEdgeColor = 'none';
h(2).MarkerFaceColor = [0 0 0]+.5;
view([45 30])
zlim([0 1])
ylim([8 32]);
%}


for val_loop=1:100
Signal = [];State=[];
Nature = [];
dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(val_loop,'%04d')]);

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

LayerNumber = 1;
numHiddenUnits = 25;

for TrainLength = 5
for GapLength = 00
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
        miniBatchSize = 16*8;
        
        RdN = 1:numel(X_RAW);%Infomation(il).randomY;
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
        acc = sum(YPred == YTest)./numel(YTest);
        
        for sl=1:numel(RdN)
            XTrain{sl,1} = X_RAW{RdN(sl),1}(:,input_Time(randperm(TrainLength)));
            YTrain(sl,1) = Y_RAW{RdN(sl),1}(1);
        end
        YTest = categorical(YTrain);
        YPred = classify(net,XTrain, ...
                'MiniBatchSize',miniBatchSize, ...
                'SequenceLength','longest');
        acc_shf = sum(YPred == YTest)./numel(YTest);
        
        
        ACCall(val_loop,TrainLength-4,GapLength+1,il) = acc;
        ACCall_shf(val_loop,TrainLength-4,GapLength+1,il) = acc_shf;
    end
end
end
end


F = figure;
subplot(1,3,1);
boxplot(ACCall);
ylim([0.7 1.0])
subplot(1,3,2);
plot([ACCall,ACCall_shf]','k');
xlim([0.5 2.5]);
ylim([0.7 1.0])
subplot(1,3,3);
boxplot(ACCall_shf);
ylim([0.7 1.0])



