clear;clc;
close all
WholePath = 'D:\ECoG\timeseries_roi_213\';
Type = {'11';'22';'33';'12';'21';'23';'31'};


% save(fullfile(WholePath,'Raw.mat'),'*_RAW');


%% pars
for LayerNumber = 1
for numHiddenUnits = 25
numHiddenUnits

for GapLength = 0

for TrainLength = 5


    TrainLength
    GapLength



    maxEpochs = 10;
    predicted_Time = 1;%62;
    input_Time = 61+(-TrainLength+1:0)-GapLength;

    for val_loop = 1:500
        %% load data
        Infomation=[];
        Signal = [];State=[];
        Nature = [];Ns0={};
        dest = fullfile(WholePath,'PCA_TEST',['Shuffle_',num2str(val_loop,'%04d')]);



        for tl = numel(Type):-1:1
            x=load(fullfile(dest,['WO_Signal',Type{tl},'.mat']));
            Y = x.Y;
            Signal = cat(3,Signal,Y);
            x=load(fullfile(dest,['WO_State',Type{tl},'.mat']));
            S = x.S;
            State = cat(3,State,S);

            N = tl*ones(size(S,3),1);
            Ns0{tl,1} = randperm(numel(N))+numel(Nature);
            Nature = cat(1,Nature,N);

        end
        X_RAW = {};
        Y_RAW = {};
        N = size(Signal,3);
        for sl=1:N
            X_RAW{sl,1} = Signal(:,:,sl)';
            Y_RAW{sl,1} = Nature(sl);%State(:,:,sl)';
        end

        Signal=[];
        x=[];



        Nr = randperm(N);

        close all
        Ns = Ns0;
        for ix=1:numel(Ns0)
            RdX = randperm(numel(Ns{ix}));
            Ns{ix} = Ns{ix}(RdX);
        end

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
        Y = Nr(1+floor(N/10*9):end);



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
        miniBatchSize = 16*8;

        layers = [];
        switch LayerNumber
            case 1
                layers = [ ...
                    sequenceInputLayer(inputSize)
                    lstmLayer(numHiddenUnits,'OutputMode','last')
                    fullyConnectedLayer(numClasses)
                    softmaxLayer
                    classificationLayer];
            case 2
                layers = [ ...
                    sequenceInputLayer(inputSize)
                    lstmLayer(numHiddenUnits,'OutputMode','last')
                    lstmLayer(numHiddenUnits,'OutputMode','last')
                    fullyConnectedLayer(numClasses)
                    softmaxLayer
                    classificationLayer];
            case 3
                layers = [ ...
                    sequenceInputLayer(inputSize)
                    lstmLayer(numHiddenUnits,'OutputMode','last')
                    lstmLayer(numHiddenUnits,'OutputMode','last')
                    lstmLayer(numHiddenUnits,'OutputMode','last')
                    fullyConnectedLayer(numClasses)
                    softmaxLayer
                    classificationLayer];
        end


        options = trainingOptions('adam', ...
            'ExecutionEnvironment','gpu', ...
            'GradientDecayFactor',0.5,...
            'GradientThresholdMethod','global-l2norm',...
            'GradientThreshold',.001, ...
            'MaxEpochs',maxEpochs, ...
            'MiniBatchSize',miniBatchSize, ...
            'ValidationData',{XValidation,YValidation},...
            'ValidationFrequency',1,...
            'Shuffle','every-epoch', ...
            'InitialLearnRate',0.01, ...
            'LearnRateSchedule','piecewise', ...
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


        Infomation(1 ).TrainingLoss = info.TrainingLoss;
        Infomation(1 ).TrainingAccuracy = info.TrainingAccuracy;
        Infomation(1 ).ValidationLoss = info.ValidationLoss;
        Infomation(1 ).TestACC = acc;
        Infomation(1 ).YTest = YTest;
        Infomation(1 ).YPred = YPred;
        Infomation(1 ).NatureTest = NatureTest;
        Infomation(1 ).net = net;
        Infomation(1 ).randomY = Y;

        filename = fullfile(dest,...
            ['InputShuffle_Infomation0_GapLength_',num2str(GapLength,'%02d'),...
            '_TrainLength_',num2str(TrainLength,'%02d'),...
            '_Layers_',num2str(LayerNumber,'%02d'),...
            '_numHiddenUnits_',num2str(numHiddenUnits,'%03d'),'.mat']);
        %save(filename,'Infomation');
        parsave(filename,Infomation)
    end

    clear *RAW X* Y*
    %
    % dest = fullfile(WholePath,'Cross-Validation-V025_3F');
    % if ~exist(dest,'dir');mkdir(dest);end



end
end
end
end




function parsave(filename,A)
save(filename,'A');
end
% for GapLength = 0:15;
% for TrainLength = 5:2:20;
% for numHiddenUnits = 50:50:300;
