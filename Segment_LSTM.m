clear;clc;
close all
%% load data
WholePath = 'H:\ECoG\timeseries_roi_213\';
Type = {'11';'22';'33';'12';'21';'23';'31'};
Signal = [];State=[];
Nature = [];Ns={};
for tl = numel(Type):-1:1
    dest = fullfile(WholePath,'ICA_TEST',['Shuffle_',num2str(1,'%04d')]);

    load(fullfile(dest,['WO_Signal',Type{tl},'.mat']));
        
    Signal = cat(3,Signal,Y);
    load(fullfile(dest,['WO_State',Type{tl},'.mat']));
    
    S=S/S(1,end,1)*tl;
    
    State = cat(3,State,S);
    
    N = tl*ones(size(S,3),1);
    Ns{tl,1} = randperm(numel(N))+numel(Nature);
    Nature = cat(1,Nature,N);
    
end
X_RAW = {};
Y_RAW = {};
N = size(Signal,3);
for sl=1:N
    X_RAW{sl,1} = Signal(:,:,sl)';
    Y_RAW{sl,1} = State(:,:,sl)';
end
clear Signal Y

%% pars

    
for numHiddenUnits = 50%[125 150 175 200] %25
    
for TrainLength = 5%3:16
    
    
    
for GapLength = 0:15
    Infomation=[]

    numHiddenUnits
    TrainLength
    GapLength
    
    %Nr = randperm(N);

    maxEpochs = 200;
    predicted_Time = 62;
    input_Time = 61+(-TrainLength+1:0)-GapLength;
    
    
    
for val_loop = 1:10
    close all
    
    X=[];
    for ix=1:numel(Ns)
        R = Ns{ix};
        X=cat(2,X,R(1:floor(numel(R)/10*8)));
    end
    V=[];
    for ix=1:numel(Ns)
        R = Ns{ix};
        V=cat(2,V,R(1+floor(numel(R)/10*8):floor(numel(R)/10*9)));
    end
    Y=[];
    for ix=1:numel(Ns)
        R = Ns{ix};
        Y=cat(2,Y,R(1+floor(numel(R)/10*9):end));
    end
    %Y = Nr(1+floor(N/10*9):end);
    
    %XV = Nr(1:floor(N/10*9));
    %Y = Nr(1+floor(N/10*9):end);    
    %Bs = floor(numel(XV)/10);
    %X = XV; X((val_loop-1)*Bs+(1:Bs)) = [];
    %V = XV((val_loop-1)*Bs+(1:Bs));
    %X = XV(1:floor(numel(XV)/10*9));
    %V = XV(1+floor(numel(XV)/10*9):end);
    
    
    %% Train data & Validation data
    XTrain = {};
    YTrain = [];
    for sl=1:numel(X)
        XTrain{sl,1} = X_RAW{X(sl),1}(:,input_Time);
        YTrain(sl,1) = Y_RAW{X(sl),1}(predicted_Time);
    end
    YTrain = categorical(YTrain);
    
    XValidation = {};
    YValidation = [];
    for sl=1:numel(V)
        XValidation{sl,1} = X_RAW{V(sl),1}(:,input_Time);
        YValidation(sl,1) = Y_RAW{V(sl),1}(predicted_Time);
    end
    YValidation = categorical(YValidation);
    
    %% LSTM
    inputSize = 100;
    numClasses = 7;
    miniBatchSize = 16;
    
    layers = [ ...
        sequenceInputLayer(inputSize)
        lstmLayer(numHiddenUnits,'OutputMode','last')
        dropoutLayer(0.25)
        fullyConnectedLayer(numClasses)
        softmaxLayer
        classificationLayer];
    
    
    options = trainingOptions('sgdm', ...
        'ExecutionEnvironment','gpu', ...
        'GradientThreshold',.1, ...
        'GradientThresholdMethod','global-l2norm',...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'ValidationData',{XValidation,YValidation},...
        'ValidationFrequency',1,...
        'SequenceLength','longest', ...
        'Shuffle','every-epoch', ...
        'InitialLearnRate',0.005, ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropPeriod',100, ...
        'LearnRateDropFactor',0.2, ...
        'Verbose',0, ...
        'Plots','training-progress');% %none
    
    [net,info] = trainNetwork(XTrain,YTrain,layers,options);
    
    %% Test
    XTest = {};
    YTest = [];
    for sl=1:numel(Y)
        XTest{sl,1} = X_RAW{Y(sl),1}(:,input_Time);
        YTest(sl,1) = Y_RAW{Y(sl),1}(predicted_Time);
    end
    YTest = categorical(YTest);
    NatureTest = Nature(Y);
    
   
    
    YPred = classify(net,XTest, ...
        'MiniBatchSize',miniBatchSize, ...
        'SequenceLength','longest');
    
    acc = sum(YPred == YTest)./numel(YTest);
    
    
    Infomation(val_loop ).TrainingLoss = info.TrainingLoss;
    Infomation(val_loop ).TrainingAccuracy = info.TrainingAccuracy;
    Infomation(val_loop ).ValidationLoss = info.ValidationLoss;
    Infomation(val_loop ).ValidationAccuracy = info.ValidationAccuracy;
    Infomation(val_loop ).TestACC = acc;
    Infomation(val_loop ).YTest = YTest;
    Infomation(val_loop ).YPred = YPred;
    Infomation(val_loop ).NatureTest = NatureTest;
    Infomation(val_loop ).net = net;
end

    dest = fullfile(WholePath,'Cross-Validation-03');
    if ~exist(dest,'dir');mkdir(dest);end
    
    filename = fullfile(WholePath,'Cross-Validation-03',...
        ['Infomation_GapLength_',num2str(GapLength,'%02d'),...
        '_TrainLength_',num2str(TrainLength,'%02d'),...
        '_numHiddenUnits_',num2str(numHiddenUnits,'%03d'),'.mat']);
    %save(filename,'Infomation');
    parsave(filename,Infomation)
end
end
end




function parsave(filename,A)
    save(filename,'A');
end
% for GapLength = 0:15;
% for TrainLength = 5:2:20;
% for numHiddenUnits = 50:50:300;
