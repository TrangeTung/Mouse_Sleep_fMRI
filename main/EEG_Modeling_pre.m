clc;clear
close all
codepath = 'D:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'D:\ECoG\power_yyl\';

Signal = [];
State = [];
for idx=1:46
    
    filepath = fullfile(WholePath,['m',num2str(idx)]);
    Sdir = dir(fullfile(filepath,'x_S5.mat'));
    if numel(Sdir)==0
        Sdir = dir(fullfile(filepath,'ch5_S*.mat'));
        if numel(Sdir)==0
            Sdir = dir(fullfile(filepath,'ch4_S*.mat'));
        end
    end
    
    if numel(Sdir)~=1
       idx
       break;
    end
    
    load(fullfile(filepath,Sdir.name));
    
    try 
        ff = f1(1:117);
        E = S1(1:14400,1:117);
    catch
        ff = f(1:117);
        E = S(1:14400,1:117);
    end
    
    
    
    load(fullfile('D:\ECoG\state_check\',['s',num2str(idx),'.mat']));
    all_state(all_state==3)=2;
    all_state(all_state==4)=3;

    Signal = [Signal;E];
    State = [State;all_state];

end

Cs = 60; % cut sec
s12_=0;s23_=0;s31_=0;s21_=0;
s11_=0;s22_=0;s33_=0;
for sl=Cs+1:numel(State)-Cs
    if State(sl)==1 & State(sl+1)==2
        if numel(unique(State(sl+(-Cs:-1))))==1
            s12_=s12_+1;XX12_(s12_)=sl;
        end
    end
    if State(sl)==2 & State(sl+1)==3
            s23_=s23_+1;XX23_(s23_)=sl;
    end
    if State(sl)==3 & State(sl+1)==1
            s31_=s31_+1;XX31_(s31_)=sl;
    end
    if State(sl)==2 & State(sl+1)==1
        if numel(unique(State(sl+(-Cs:-1))))==1
            s21_=s21_+1;XX21_(s21_)=sl;
        end
    end
    
    y = State(sl+(-Cs:3));
    y_ = unique(y);
    if numel(y_)==1;
        if y_==1; s11_=s11_+1;XX11_(s11_)=sl;end
        if y_==2; s22_=s22_+1;XX22_(s22_)=sl;end
        if y_==3; s33_=s33_+1;XX33_(s33_)=sl;end
    end
    
end


clear Y*_ S 
Type = {'11';'22';'33';'12';'21';'23';'31'};
for tl=1:numel(Type)
    bins = 1;
    if tl==3 | tl==7| tl==6;bins=1;end
    Y=[];S=[];
    eval(['sx=s',Type{tl},'_;']);
    eval(['XX=XX',Type{tl},'_;']);
    
    BINS = [165;165;180;170;165;170;165];
    
    if tl==1;StateIndx = 1*ones(Cs*2+1,1);end
    if tl==2;StateIndx = 2*ones(Cs*2+1,1);end
    if tl==3;StateIndx = 3*ones(Cs*2+1,1);end
    
    if tl==4;StateIndx = [1*ones(Cs+1,1);2*ones(Cs,1)];end
    if tl==5;StateIndx = [2*ones(Cs+1,1);1*ones(Cs,1)];end
    if tl==6;StateIndx = [2*ones(Cs+1,1);3*ones(Cs,1)];end
    if tl==7;StateIndx = [3*ones(Cs+1,1);1*ones(Cs,1)];end
    
    
    for loop=1:BINS(tl)
        M=[]; 
        XX = XX(randperm(numel(XX)));
        for sl=1:bins
            M(:,:,sl) = Signal(XX(sl)+(-Cs:Cs),:);
        end
       Y(:,:,loop) = mean(M,3);
       
       S(1,:,loop) = StateIndx;
    end
    
    eval(['Y',Type{tl},'_=Y;']);
    
%     filename = fullfile(WholePath,['Signal',Type{tl},'.npy']);
%     writeNPY(Y, filename);
    save(fullfile(WholePath,['EEG_Signal',Type{tl},'.mat']),'Y');
%     filename = fullfile(WholePath,['State',Type{tl},'.npy']);
%     writeNPY(S, filename);
    save(fullfile(WholePath,['EEG_State',Type{tl},'.mat']),'S');
end


1;





