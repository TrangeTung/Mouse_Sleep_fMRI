function F = MY_get_figure_Displacement_Maximum(rp,legend_txt,title_txt)
%% rp : N¡Á6 (x,y,z,x',y',z')
% F: Figure
F = figure;
suptitle(title_txt);
set(gca,'FontWeight','bold');
set(F,'position',[100 100 550 250]);

color = {[138 043 226];[002 152 182];[135 199 000];[248 199 000];[005 139 085];[000 000 255]};

%% rp
for count = 1:numel(legend_txt)
    rp_Var = rp(:,count);
    for k = 1:numel(rp_Var)/5
        rp_subVar = rp_Var((k-1)*5+1:(k-1)*5+5);
        plot((count*3-1+k*0.1)*ones(numel(rp_subVar),1),rp_subVar,'o','color',color{2}/255,'MarkerFaceColor',color{2}/255);
        hold on;
    end
    
    rp_Var_mean = mean(rp_Var);
    rp_Var_SEM = std(rp_Var);%/sqrt(numel(rp_Var));
    count_center = count*3-1+numel(rp_Var)/10*0.1;
    count_label(count) = count_center;
    plot(count_center-9*0.05:0.05:count_center+9*0.05,rp_Var_mean*ones(19,1),'k','LineWidth',2)
    plot((count_center-8*0.05)*ones(1/0.05+1,1),rp_Var_mean*(0:0.05:1),'k','LineWidth',2)
    plot((count_center+8*0.05)*ones(1/0.05+1,1),rp_Var_mean*(0:0.05:1),'k','LineWidth',2)
    plot(count_center-3*0.05:0.05:count_center+3*0.05,(rp_Var_mean+rp_Var_SEM)*ones(7,1),'k','LineWidth',2)
    plot(count_center-3*0.05:0.05:count_center+3*0.05,(rp_Var_mean-rp_Var_SEM)*ones(7,1),'k','LineWidth',2)
    step = rp_Var_SEM*2/10;
    plot(count_center*ones(11,1),(rp_Var_mean-rp_Var_SEM):step:(rp_Var_mean+rp_Var_SEM),'k','LineWidth',2)
end
hold off;
if contains(title_txt,'Translation')
    ylabel('mm');
    ymax = 0.04;
else
    ylabel('degree');
    ymax = 2;
end
set(gca,'FontSize',10,'FontWeight','bold');
set(gca, 'LineWidth',1.5);
box off;

axis([0 count*3+3 0 ymax])
xtextp = count_label;
ytextp = -ymax/25*ones(numel(xtextp),1);
set(gca,'xticklabel','');
text(xtextp,ytextp,legend_txt,'HorizontalAlignment','center','fontsize',10); 
xticks(xtextp)

text(max(xtextp(:)),-ymax/25*2,'(Mean¡Àstd)','FontWeight','bold')
end
