function [FD_GS_r,F] = MY_get_figure_rp_GSintensity(filepath,rp,Resp,G_Signal,mask)
%%
% rp : head motion parameters (6) (relative URL)
% G_Signal: the whole brain signal name (relative URL) {cell}
% F: Figure
color = {[138 043 226];[112 173 071];[138 043 226];[068 114 196];[237 125 049];[248 192 000];};
F = figure;
set(F,'position',[0 0 210 17*(2+3*numel(G_Signal))]*3.5);
axis normal;
%% rp
positionVector=[0.05 1-1/(2+3*numel(G_Signal)) 0.85 0.75/(2+3*numel(G_Signal))];
subplot('Position',positionVector);
cd(filepath)
rp_6 = load(rp);
rp_6(:,4:6) = rp_6(:,4:6).*pi*180;
rp_6(:,1:3) = rp_6(:,1:3);
plot(rp_6)
set(gca,'LooseInset',get(gca,'TightInset'))
xlim([1 length(rp_6)]);
ylim([min(rp_6(:))*1.1 max(rp_6(:))*1.1]);
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');box off
ylabel({'rp';'mm/ бу'},'FontSize',8,'FontWeight','bold','color','k');
xloc = length(rp_6)*1.005;
yloc = min(rp_6(:))*1.1*0.15+max(rp_6(:))*1.1*0.85;
text(xloc,yloc,'Max/Min','HorizontalAlignment','left','fontsize',7,'color',color{1}/255);
yloc = min(rp_6(:))*1.1*0.35+max(rp_6(:))*1.1*0.65;
text(xloc,yloc,[num2str(max(max(rp_6(:,1:3))),'%0.2f'),'/',num2str(min(min(rp_6(:,1:3))),'%0.2f')],...
    'HorizontalAlignment','left','fontsize',7,'color',color{1}/255);
yloc = min(rp_6(:))*1.1*0.55+max(rp_6(:))*1.1*0.45;
text(xloc,yloc,[num2str(max(max(rp_6(:,4:6))),'%0.2f'),'/',num2str(min(min(rp_6(:,4:6))),'%0.2f')],...
    'HorizontalAlignment','left','fontsize',7,'color',color{1}/255);
%% Resp
% fs = 100;
% [s_0]=textread(Resp{1},'%s');
% s_0(1:81)=[];NUM_ram = str2double(s_0);
% NUM = reshape(NUM_ram(:),[8 numel(NUM_ram)/8])';
% % find trigger
% trigger = NUM(:,4);
% trigger_norm = imbinarize(trigger,(max(trigger)+min(trigger))/2);
% RAM = find(diff(trigger_norm)<0);
% % resp
% resp = NUM(RAM(3):end,6);
% [b,a] = butter(10,0.3);
% resp_ram = filtfilt(b,a,resp);
% [~,LOCS] = findpeaks(resp_ram);
% count = [1:Resp{2}*fs:numel(resp_ram)+Resp{2}*fs*2];
% for ram = 1:length(rp_6)
%     Resp_rate(ram) = numel(find(LOCS<count(ram+2)&LOCS>count(ram)))*(60/(Resp{2}*2));
% end
% if Resp{2}<1;[b,a] = butter(10,(round(Resp{2}/1.2*10))/10);Resp_rate = filtfilt(b,a,Resp_rate);end

%% FD
positionVector=[0.05 1-2/(2+3*numel(G_Signal)) 0.85 1/(2+3*numel(G_Signal))];
subplot('Position',positionVector);
cd(filepath)

rp_6 = load(rp);
rp_6 = rp_6 - [zeros(1,6);rp_6(1:end-1,:)];
[~,rp_6] = gradient(rp_6);
rp_6(:,4:6) = rp_6(:,4:6).*pi./180.*5; %% mm
rp_6(:,1:3) = rp_6(:,1:3)./6;
rp_displace = sqrt(rp_6(:,1).^2+rp_6(:,2).^2+rp_6(:,3).^2+rp_6(:,4).^2+rp_6(:,5).^2+rp_6(:,6).^2);

plot(rp_displace,'color',color{5}/255)
set(gca,'LooseInset',get(gca,'TightInset'))
xlim([1 length(rp_6)]);
ylim([0 max(rp_displace(:))*1.1]);
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');box off;
ylabel({'FD';'\mum'},'FontSize',8,'FontWeight','bold','color','k');
xloc = length(rp_6)*1.005;
yloc = min(rp_displace(:))*1.1*0.00+max(rp_displace(:))*1.1*1.00;
text(xloc,yloc,'FD','HorizontalAlignment','left','fontsize',7,'color',color{5}/255);
yloc = min(rp_displace(:))*1.1*0.15+max(rp_displace(:))*1.1*0.85;
text(xloc,yloc,[num2str(max(rp_displace(:)*1000),'%0.2f'),'/',num2str(min(rp_displace(:)*1000),'%0.2f')],...
    'HorizontalAlignment','left','fontsize',7,'color',color{1}/255);
% hold on;plot((Resp_rate-min(Resp_rate)*0.9)/max((Resp_rate-min(Resp_rate)*0.9))*max(rp_displace(:)),'color',color{4}/255);hold off;
% yloc = min(rp_displace(:))*1.1*0.30+max(rp_displace(:))*1.1*0.70;
% text(xloc,yloc,'Resp rate','HorizontalAlignment','left','fontsize',7,'color',color{4}/255);
% yloc = min(rp_displace(:))*1.1*0.45+max(rp_displace(:))*1.1*0.55;
% text(xloc,yloc,[num2str(round(max(Resp_rate(:)))),'/',num2str(round(min(Resp_rate(:))))],...
%     'HorizontalAlignment','left','fontsize',7,'color',color{1}/255);
% yloc = min(rp_displace(:))*1.1*0.75+max(rp_displace(:))*1.1*0.25;
% text(xloc,yloc,'FD&Resp','HorizontalAlignment','left','fontsize',7,'color','k');
% yloc = min(rp_displace(:))*1.1*0.95+max(rp_displace(:))*1.1*0.05;
% text(xloc,yloc,['r=',num2str(corr(rp_displace(:),Resp_rate(:)),'%0.2f')],'HorizontalAlignment','left','fontsize',7,'color','k');

for kk = 1:numel(G_Signal)
    nii = load_nii(G_Signal{kk});
    GS = fmask(nii.img,mask);
    mask_RAM = mean(GS,2);
    GS(mask_RAM<=0,:) = [];
    GS = double(GS)-mean(GS,2);
    GS(isnan(GS)) = 0;
    
    positionVector=[0.05 1-(3*kk-0.25)/(2+3*numel(G_Signal)) 0.85 0.75/(2+3*numel(G_Signal))];
    subplot('Position',positionVector);
    plot(mean(GS,1),'k')
    xlim([1 numel(mean(GS,1))]);
    ylim([min(mean(GS,1))*1.1 max(mean(GS,1))*1.1]);
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');box off
    ylabel({'GS';'A.U.'},'FontSize',8,'FontWeight','bold','color','k');
    set(gca,'LooseInset',get(gca,'TightInset'))
    yloc = min(mean(GS,1))*1.1*0.25+max(mean(GS,1))*1.1*0.75;
    text(xloc,yloc,[num2str(max(mean(GS,1)),'%0.2f'),'/',num2str(min(mean(GS,1)),'%0.2f')],...
        'HorizontalAlignment','left','fontsize',7,'color',color{1}/255);
    yloc = min(mean(GS,1))*1.1*0.55+max(mean(GS,1))*1.1*0.45;
    text(xloc,yloc,'FD&GS','HorizontalAlignment','left','fontsize',7);
    yloc = min(mean(GS,1))*1.1*0.85+max(mean(GS,1))*1.1*0.15;
    x = rp_displace(:);
    y = mean(GS,1);
    FD_GS_r(kk) = corr(x(5:end),y(5:end)');
    text(xloc,yloc,['r=',num2str(FD_GS_r(kk),'%0.4f')],'HorizontalAlignment','left','fontsize',7);
    positionVector=[0.05 1-(3*kk+2)/(2+3*numel(G_Signal)) 0.85 2.25/(2+3*numel(G_Signal))];
    subplot('Position',positionVector);
    imagesc(GS);box off;
    set(gca,'ytick',[],'xtick',[],'xcolor','w','ycolor','w')
    ylabel({G_Signal{kk};[num2str(size(GS,1)),' voxels']},'FontSize',8,'FontWeight','bold','Interpreter','none','color','k');
    caxis([mean(GS(:))-2*std(GS(:)) mean(GS(:))+2*std(GS(:))]);
    x = [size(GS,2)*0.05 size(GS,2)*0.05 size(GS,2)*0.05+10 size(GS,2)*0.05+10];
    y = [size(GS,1)*0.1 size(GS,1)*0.1*2 size(GS,1)*0.1*2 size(GS,1)*0.1];
    patch(x,y,'k','FaceAlpha',.75,'EdgeColor','none')
    text(size(GS,2)*0.05+11,size(GS,1)*0.1+size(GS,1)*0.1/2,'10 TR','HorizontalAlignment','left','Fontweight','bold','color','k')
end

end