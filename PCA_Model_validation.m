%% Power validation
close all
rl=0; sl=0; nl=0;
Ac=zeros(61,2393,200);
Rc=zeros(61,2393,200);
Nc=zeros(61,2393,200);
for idx = 1:19
    
    load(fullfile('H:\ECoG\lfp1_500hz',['x_S_500_',num2str(idx),'.mat']));
    power = 10*log10(S1);
    %power = medfilt2(power,[1 50]);
    power(:,f1<52 & f1>48)=nan;
    power = fillmissing(power','linear')';
    
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx+27),'.mat']));
    
    X = all_state(5:end);
    dX = [X-[X(1);X(1:end-1)]];
    Cut=30;
    
    Event = power;
     %figure;plot(f1,mean(power,1))
    locs = find(dX==-2); locs(locs<Cut|locs>numel(X)-Cut)=[];
    for loop=1:numel(locs); rl=rl+1;
        Ac(:,:,rl) = Event(locs(loop)+(-Cut:Cut),:);
    end
    locs = find(dX==1); locs(locs<Cut|locs>numel(X)-Cut)=[];
    for loop=1:numel(locs); sl=sl+1;
        Rc(:,:,sl) = Event(locs(loop)+(-Cut:Cut),:);
    end
    
    y = (dX==0 & X==3); 
    SE = strel('disk',3);
    y_ = imerode(y,SE); 
    locs = find(y_==1); locs(locs<Cut|locs>numel(X)-Cut-6)=[];
    LL=locs; locs = LL(1:Cut:end);
    for loop=1:numel(locs); nl=nl+1;
        Nc(:,:,nl) = Event(locs(loop)+(-Cut:Cut),:);
    end
end
Ac(:,:,rl+1:end)=[];
Rc(:,:,sl+1:end)=[];
Nc(:,:,nl+1:end)=[];



t=-30:30;
F = figure;
subplot(1,2,1)
X=mean(Ac,3)-mean(Nc,3);
x_ = reshape(Ac,[],size(Ac,3));
y_ = reshape(Nc,[],size(Nc,3));
[~,p] = ttest2(x_',y_','tail','right');
pval = reshape(p,size(X));
X(pval>0.05)=nan;
imagesc(t,f1,X');axis xy;
colorbar; title('Spectrogram');
colormap('jet');caxis([-2 2])
subplot(1,2,2)
X=mean(Rc,3)-mean(Nc,3);
x_ = reshape(Rc,[],size(Rc,3));
y_ = reshape(Nc,[],size(Nc,3));
[~,p] = ttest2(x_',y_','tail','right');
pval = reshape(p,size(X));
X(pval>0.15)=nan;
imagesc(t,f1,X');axis xy;
colorbar; title('Spectrogram');
colormap('jet');caxis([-5 5])


%% Event Validation MRI scanner
%{
rl=0; sl=0; nl=0;
Ac=[]; Rc=[]; Nc=[];
for idx = [28:34 37:39 41:46]
    
    load(fullfile('H:\ECoG\event_test','ripple_peak',['ripple_peak',num2str(idx),'.mat']));
    x = nrem_peak; RIPPLE = mean(reshape(x,1024,[]),1)*100; 
    RIPPLE = envelope(RIPPLE,7,'rms'); RIPPLE=RIPPLE-mean(RIPPLE(1:10)); RIPPLE(RIPPLE<0)=0;
    %RIPPLE(RIPPLE>50/1024)=1;
    
    load(fullfile('H:\ECoG\event_test','spindle',['spindle_peak',num2str(idx),'.mat']));
    x = spindle; SPINDLE = mean(reshape(x,1024,[]),1)*100; 
    SPINDLE = envelope(SPINDLE,7,'rms'); SPINDLE=SPINDLE-mean(SPINDLE(1:10)); SPINDLE(SPINDLE<0)=0;
    %SPINDLE(SPINDLE>50/1024)=1; SPINDLE(SPINDLE~=1)=0;
    
    load(fullfile('H:\ECoG\state_check\',['s',num2str(idx),'.mat']));
    
    X = all_state;
    dX = [X-[X(1);X(1:end-1)]];
    Cut=30;
    
    Event = SPINDLE;
    
    locs = find(dX==-2); locs(locs<Cut|locs>numel(X)-Cut)=[];
    for loop=1:numel(locs); rl=rl+1;
        Ac(:,rl) = Event(locs(loop)+(-Cut:Cut));
    end
    locs = find(dX==1); locs(locs<Cut|locs>numel(X)-Cut)=[];
    for loop=1:numel(locs); sl=sl+1;
        Rc(:,sl) = Event(locs(loop)+(-Cut:Cut));
    end
    
    y = (dX==0 & X==3); 
    SE = strel('disk',Cut);
    y_ = imerode(y,SE); 
    locs = find(y_==1); locs(locs<Cut|locs>numel(X)-Cut)=[];
    LL=locs; locs = LL(1:Cut:end);
    for loop=1:numel(locs); nl=nl+1;
        Nc(:,nl) = Event(locs(loop)+(-Cut:Cut));
    end
end

x = mean(Ac,2); shfs=min(find(x==min(x)))-Cut-1; Ac=[Ac(shfs+1:end,:);Ac(end,:).*ones(shfs,size(Ac,2))];
x = mean(Rc,2); shfs=min(find(x==min(x)))-Cut-1; Rc=[Rc(shfs+1:end,:);Rc(end,:).*ones(shfs,size(Rc,2))];

close all
F = figure('Position', [680 541 451 437]);
x = -Cut:Cut;  
Ms = mean(Ac,2); SEMs = std(Ac,0,2)/sqrt(size(Ac,2));
plot(x,Ms,'k'); hold on;
patch('XData',[x,fliplr(x)],'YData',[Ms+SEMs;flipud(Ms-SEMs)],...
    'EdgeColor','none','FaceColor','y','FaceAlpha',0.2);
Ms = mean(Rc,2); SEMs = std(Rc,0,2)/sqrt(size(Rc,2));
plot(x,Ms,'k')
patch('XData',[x,fliplr(x)],'YData',[Ms+SEMs;flipud(Ms-SEMs)],...
    'EdgeColor','none','FaceColor','b','FaceAlpha',0.2);
Ms = mean(Nc,2); SEMs = std(Nc,0,2)/sqrt(size(Nc,2));
plot(x,Ms,'k')
patch('XData',[x,fliplr(x)],'YData',[Ms+SEMs;flipud(Ms-SEMs)],...
    'EdgeColor','none','FaceColor','G','FaceAlpha',0.2);
yl=ylim;ylim([0 max(yl)]);
xlim([-Cut Cut]);
xlabel('Time (s)');
ylabel('Probability (%)');
set(gca,'box','off','tickdir','out','fontsize',15,...
    'linewidth',1.5,'ticklength',[0.02 0.1])

clear p*
for ix=1:61
   x0 = Nc(ix,:);
   y1 = Ac(ix,:);
   y2 = Rc(ix,:);
   [~,p1(ix)] = ttest2(x0(:),y1(:),'tail','right');
   [~,p2(ix)] = ttest2(x0(:),y2(:),'tail','right');
end

%}
