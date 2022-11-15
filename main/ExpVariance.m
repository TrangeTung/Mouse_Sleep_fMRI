WholePath = 'F:\ECoG\pc_yyl\';

R2all = cell(1,46);
for idx=[1:41 43:46]
    cd(fullfile(WholePath,num2str(idx,'%02d')));
    coeff = load('coeff.txt');
    score = load('PC_all.txt');
    
    ySC = score*coeff';
    R2 = zeros(100,1);
    for vl = 1:100
        DMx = score(:,vl)*coeff(:,vl)';
        Vec = DMx(:)'; VecP = ySC(:)';
        mu = mean(VecP,2).*mean(Vec,2) - mean(Vec.*VecP,2);
        md = mean(VecP,2).^2 - mean(VecP.^2,2);
        b = mean(Vec,2)-mu./md.*mean(VecP,2);
        Vecf = mu./md.*VecP + b;
        SSR = sum((Vecf-mean(Vec,2)).^2,2);
        SST = sum((Vec-mean(Vec,2)).^2,2);
        R2(vl) = SSR./SST;
    end
    
    R2all{idx} = cumsum(R2/sum(R2));
end
R2all(42)=[];


F = figure;
V = 100*cell2mat(R2all(1:27))';
Ms = mean(V,1); x = 1:size(V,2);
SEM = std(V,0,1)/sqrt(size(V,1));
semilogx(x,Ms,'b','linewidth',.5);
hold on;patch([x,fliplr(x)],[Ms-SEM,fliplr(Ms+SEM)],...
    ones(size([Ms-SEM,fliplr(Ms+SEM)])),'edgecolor',...
    'none','facecolor','b','facealpha',.3);
V = 100*cell2mat(R2all(28:end))';
Ms = mean(V,1);
SEM = std(V,0,1)/(size(V,1));
plot(x,Ms,'r','linewidth',.5);
hold on;patch([x,fliplr(x)],[Ms-SEM,fliplr(Ms+SEM)],...
    ones(size([Ms-SEM,fliplr(Ms+SEM)])),'edgecolor',...
    'none','facecolor','r','facealpha',.3);
set(gca,'fontsize',15,'linewidth',.5,'box','off','tickdir','out');
ylim([0 100]);


