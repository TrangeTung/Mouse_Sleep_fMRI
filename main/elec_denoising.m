%% load data
cd(path);
ecgAll=load('ch1.mat');
fsEcg=load('fsEcg.mat');
wavAll=load('trigger.mat');
fsWav=load('Fs.mat');
x=1/fsEcg:1/fsEcg:size(ecgAll,1)/fsEcg;
figure;
plot(x,ecgAll(:,1)); 
%% 
secStart=85;%change
save('secStart.mat','secStart');
numTriggerDum=22*4;%change
for c=1:1:2 % determine true secStart, after dummy scan. 
    %10 seconds before true EPI(dummy scan excluded)
%ecg=ecgAll(secStart*fsEcg:end,:);
wav=wavAll(secStart*fsWav:end,:);
trigger=wav;
difT=diff(trigger);%≤Ó÷µ
startTrigger=find(difT<-0.5);
invalidTrigger=[];
a=2;
b=1;
numTrigger=size(startTrigger,1);
while(a<=size(startTrigger,1))
    if isempty(find(difT(startTrigger(a-1):startTrigger(a),1)>0.5, 1))
        invalidTrigger(b,1)=a;
        invalidTrigger(b,2)=startTrigger(a);
        invalidTrigger(b,3)=startTrigger(a-1);
        if difT(startTrigger(a))>difT(startTrigger(a-1))
        startTrigger(a)=[];
        else
            startTrigger(a-1)=[];
        end
        b=b+1;
    end
    a=a+1;
end
% if c==1
% %secStart=secStart+startTrigger(numTriggerDum+1)/fsWav-10;%8 seconds before 
% secStart=secStart-10;
%true EPI(dummy scan excluded)
%end
end
secStart0=secStart;
save('secStart0.mat','secStart0');
save('startTrigger.mat','startTrigger');
clear ecgAll;
%% para
ecg=ecgAll(3,secStart0*fsEcg:end);
EEG.srate=fix(data.streams.EOG3.fs);
EEG.pnts=length(ecg);
nSlice=22;
i=1;
for a=1:nSlice:size(startTrigger,1)
    trig(i)=startTrigger(a,1)/fsWav*fsEcg;
    i=i+1;
end
trig=ceil(trig);
clear data;
%% denoising
for i=[1:16]
    EEG.data=ecg(i,:);
    de_data=fmrib_fastr(EEG,0,1,30,trig,0,0,0);
    name_res=['ch',num2str(i+2),'.mat'];
    save(name_res,'de_data','-v7.3');
    clear de_data;
end