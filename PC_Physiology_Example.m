clc;clear;

% example
%{
WholePath = 'F:\ECoG\m10_2\';
cd(WholePath)
hdr = spm_vol('2dseq.nii');
img = spm_read_vols_4D(hdr);

lmask = spm_read_vols(spm_vol('EPI_mask.nii'));
V = fmask(img,lmask);
V = (V-mean(V,2))./std(V,0,2);
[~,score] = pca(V');

load('Resp.mat');
[b,a] = butter(3,0.005);
RX = filtfilt(b,a,Resp);
% [~,LOCS]=findpeaks(RX(1:10^5),'MinPeakDistance',Resp_fs/5);
load('secStart0.mat');
DummyTime = secStart0+8;
s = round(DummyTime*Resp_fs);

for sl = 1:7200
   sf = s + round(sl* Resp_fs);
   cs = sf + round([-20 40]*Resp_fs);
    [~,LOCS]=findpeaks(RX(min(cs):max(cs)),'MinPeakDistance',Resp_fs/5);
    RR(sl)=numel(LOCS);
end
CC=corr(score,RR');
[~,sx]=sort(CC,'descend');

score1 = (score-mean(score,2)) ./std(score,0,2);
figure;subplot(7,1,1);plot(RR);xlim([0 7200]);
subplot(7,1,2:4);imagesc(V(sx([1:2000 2001:10:end-2000 end-2000:end]),:));caxis([-2 2]);colormap('gray')
subplot(7,1,5:7);plot(-score(:,1:10)+(1:10)*-40);xlim([0 7200]);
%}







WholePath = 'F:\ECoG\pc_yyl\';
clear CC CCs
R2 = nan(46,1);
CC = nan(46,100);
for idx=1:46
    cd(fullfile(WholePath,num2str(idx,'%02d')));
    
    
    
    
    if exist('Resp.mat','file')
            PCs = load('PC_all.txt');
            
            Resp_fs = 1.017252624511719e+03;
            load('Resp.mat');
            [b,a] = butter(3,0.005);
            RX = filtfilt(b,a,Resp);
            % [~,LOCS]=findpeaks(RX(1:10^5),'MinPeakDistance',Resp_fs/5);
            load('secStart0.mat');
            DummyTime = secStart0+8;
            s = round(DummyTime*Resp_fs);

            RR = zeros(7200,1);
            for sl = 1:7200
               sf = s + round(sl* Resp_fs);
               cs = sf + round([-20 40]*Resp_fs);
                [~,LOCS]=findpeaks(RX(min(cs):max(cs)),'MinPeakDistance',Resp_fs/5);
                RR(sl)=numel(LOCS);
            end    
            save('RR.mat','RR');

            [~,~,r] = regress(RR,[ones(size(RR)),PCs(:,1:40)]);
            
            Vec = RR(:)'-r(:)'; VecP = RR(:)';
            mu = mean(VecP,2).*mean(Vec,2) - mean(Vec.*VecP,2);
            md = mean(VecP,2).^2 - mean(VecP.^2,2);
            b = mean(Vec,2)-mu./md.*mean(VecP,2);
            Vecf = mu./md.*VecP + b;
            SSR = sum((Vecf-mean(Vec,2)).^2,2);
            SST = sum((Vec-mean(Vec,2)).^2,2);
            R2(idx,1) = SSR./SST;
            CC(idx,:) = corr(PCs(:,1:100),RR);
    end

end


F = figure;
V = abs(CC(1:19,:));
Ms = mean(V,1); x = 1:size(V,2);
SEM = std(V,0,1)/sqrt(size(V,1));
semilogx(x,Ms,'b','linewidth',.5);
hold on;patch([x,fliplr(x)],[Ms-SEM,fliplr(Ms+SEM)],...
    ones(size([Ms-SEM,fliplr(Ms+SEM)])),'edgecolor',...
    'none','facecolor','b','facealpha',.3);
V =  abs(CC(28:end,:));
Ms = mean(V,1);
SEM = std(V,0,1)/(size(V,1));
plot(x,Ms,'r','linewidth',.5);
hold on;patch([x,fliplr(x)],[Ms-SEM,fliplr(Ms+SEM)],...
    ones(size([Ms-SEM,fliplr(Ms+SEM)])),'edgecolor',...
    'none','facecolor','r','facealpha',.3);
set(gca,'fontsize',15,'linewidth',.5,'box','off','tickdir','out');
ylim([0 0.4]);

F = figure;
h=histogram(R2(1:19,:),20,'BinLimits',[0 1],'Normalization','count',...
    'displaystyle','stair','EdgeColor',[0 0 1]);
hold on;
h=histogram(R2(28:end,:),20,'BinLimits',[0 1],'Normalization','count',...
    'displaystyle','stair','EdgeColor',[1 0 0]);
set(gca,'box','off','linewidth',.5,'tickdir','out','fontsize',12);







