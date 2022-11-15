function [F,post_hoc] = MY_get_figure_Displacement_Variability(rp,legend_txt,title_txt)
%% rp : cell{rp} ; Multiple Headmotion datasets
% F: Figure
% c : post-hoc results
F = figure;
suptitle(title_txt);
set(gca,'FontWeight','bold');
set(F,'position',[100 100 550 550]);

color = {[138 043 226];[002 152 182];[135 199 000];[248 199 000];[005 139 085];[000 000 255]};

%% rp
for count = 1:numel(legend_txt)
    rp_Var = rp(count,:);
    plot(count*ones(numel(rp_Var),1),rp_Var,'o','color',color{count}/255)
    hold on;
    rp_Var_mean = mean(rp_Var);
    rp_Var_SEM = std(rp_Var)/sqrt(numel(rp_Var));
    plot(count-4*0.05:0.05:count+4*0.05,rp_Var_mean*ones(9,1),'k','LineWidth',2)
    plot(count-2*0.05:0.05:count+2*0.05,(rp_Var_mean+rp_Var_SEM)*ones(5,1),'k','LineWidth',2)
    plot(count-2*0.05:0.05:count+2*0.05,(rp_Var_mean-rp_Var_SEM)*ones(5,1),'k','LineWidth',2)
    step = rp_Var_SEM*2/10;
    plot(count*ones(11,1),(rp_Var_mean-rp_Var_SEM):step:(rp_Var_mean+rp_Var_SEM),'k','LineWidth',2)
end
%% post-hoc analysis
data = rp';
[p,~,stats]=anova1(data);
close;
[c,~,h,~]=multcompare(stats,'alpha',0.05,'ctype','tukey-kramer');
close(h);
post_hoc = c;

ylabel('tSNR (A.U.)');
set(gca,'FontSize',10,'FontWeight','bold');
set(gca, 'LineWidth',1.5);
box off;
axis([0 count+1 0 10])
set(gca,'xticklabel',[{''}; legend_txt; {''}]);


% h = legend(legend_txt{:});
% set(h,'Box','off','location','northeastoutside');

end
