clc;clear
close all
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'D:\ECoG\timeseries_roi_213\';

Labels_Excel = fullfile(WholePath,'Allen_mouse_SC_213.xlsx');
[~,~,CellData_213] = xlsread(Labels_Excel);

Labels_Excel = fullfile(WholePath,'timeseries_roi_213.xlsx');
[~,~,CellData] = xlsread(Labels_Excel);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
MajorRegion = cat(1,ExpTable.abbre);


clear Idx
T = CellData_213(2:end,1);
for tl=1:numel(T)
    A = T{tl};
    B = cellfun(@(x) strcmp(x,A), MajorRegion);
    
    C = find(B==1);
    if isempty(C);
        C=133;
        tl
        A
    end
    Idx(tl)=C;
end

Ix = ExpTable.Include(Idx);
M0 = cell2mat(CellData_213(2:end,2:end));
Iall = M0(Ix==1,Ix==1);
I = -log(Iall);
Region = T(Ix==1);
I(I==-Inf)=nan; I(I>15)=nan; I = I-min(I(:));
I(isnan(I))=100;

% csvwrite(fullfile(WholePath,'Weight.csv'),I);

Signal = [];
State = [];
for idx=1:46
    %% ICA
    %
    load(fullfile(WholePath,'ICA_TEST','TimeSeries',['tc_',num2str(idx,'%02d'),'.mat']));
    beta = tc;
    beta = (beta-mean(beta,2))./std(beta,0,2);
    E = beta';
    %}
    
    
    %% PCA
    %{
    load(fullfile('H:\ECoG\pca',num2str(idx,'%02d'),'GroupPCs.mat'));
    beta = (beta-mean(beta,2))./std(beta,0,2);
    E = beta';
    %}
    
    [b,a]=butter(3,0.1);
    for ix=1:size(E,2); F(:,ix) = filtfilt(b,a,E(:,ix)'); end
    
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));
    b1 = interp1(1:14400,all_state,2:2:14400,'nearest')';
    b1(b1==3)=2;
    b1(b1==4)=3;
    
    Signal = [Signal;E];
    State = [State;b1];

end

Cs = 60; % cut TR
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

for shf_loop = 100:500
clear Y*_ S 
Type = {'11';'22';'33';'12';'21';'23';'31'};
for tl=1:numel(Type)
    bins = 1;

    Y=[];S=[];
    eval(['sx=s',Type{tl},'_;']);
    eval(['XX=XX',Type{tl},'_;']);
    
    BINS = [73;73;73;73;73;73;73];
    BINS = [165;165;180;170;165;170;165]*5;
    
    if tl==1;StateIndx = 1*ones(Cs*2+1,1);end
    if tl==2;StateIndx = 2*ones(Cs*2+1,1);end
    if tl==3;StateIndx = 3*ones(Cs*2+1,1);end
    
    if tl==4;StateIndx = [1*ones(Cs+1,1);2*ones(Cs,1)];end
    if tl==5;StateIndx = [2*ones(Cs+1,1);1*ones(Cs,1)];end
    if tl==6;StateIndx = [2*ones(Cs+1,1);3*ones(Cs,1)];end
    if tl==7;StateIndx = [3*ones(Cs+1,1);1*ones(Cs,1)];end
    
    if tl==4;StateIndx = [1*ones(Cs+1,1);4*ones(Cs,1)];end
    if tl==5;StateIndx = [2*ones(Cs+1,1);5*ones(Cs,1)];end
    if tl==6;StateIndx = [2*ones(Cs+1,1);6*ones(Cs,1)];end
    if tl==7;StateIndx = [3*ones(Cs+1,1);7*ones(Cs,1)];end
    
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
    
    dest = fullfile(WholePath,'ICA_TEST',['Shuffle_',num2str(shf_loop,'%04d')]);
    if ~exist(dest,'dir');mkdir(dest);end
    
    save(fullfile(dest,['WO_Signal',Type{tl},'.mat']),'Y');
    save(fullfile(dest,['WO_State',Type{tl},'.mat']),'S');
end
end
1;





