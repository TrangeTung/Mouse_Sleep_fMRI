function F = MY_get_figure_Frame_wise_Displacement(rp,stimulus,legend_txt)
%% rp : cell{rp} ; Multiple Headmotion datasets
% F: Figure

F = figure;
suptitle('Frame-wise Displacement');
set(gca,'FontWeight','bold');
set(F,'position',[0 0 900 900]);

color = {[138 043 226];[002 152 182];[135 199 000];[248 199 000];[005 139 085];[000 000 255]};

%% rp
for count = 1:numel(rp)
    TimeCourse_mean = rp{count}; 
    number = size(TimeCourse_mean,2);
    Time_series = mean(TimeCourse_mean,2);
    upper = Time_series+std(TimeCourse_mean,0,2)/sqrt(number);
    lower = flipud(Time_series-std(TimeCourse_mean,0,2)/sqrt(number));
    y_min = 0.2;
    y_max = 0.5;
    
    subplot(numel(rp),1,count);
    f = plot([1:numel(Time_series)],Time_series,'color',color{count}/255,'LineWidth',1);
    x = [1:numel(Time_series),fliplr(1:numel(Time_series))]';
    y = [upper;lower];
    patch(x,y,color{count}/255,'FaceAlpha',.2,'EdgeColor','none');  
    hold on;
    
    %% stimulus
    stimulus = stimulus(:);
    stimulus(stimulus~=0) = 1;
    stimulus_dev = stimulus - [stimulus(2:end);stimulus(1)];
    stimulus_up = find(stimulus_dev==1);
    stimulus_down = find(stimulus_dev==-1);
    patch_x = [stimulus_down stimulus_down stimulus_up stimulus_up]';
    patch_y = [repmat(y_min,1,numel(stimulus_down));repmat(y_max,1,numel(stimulus_up));
               repmat(y_max,1,numel(stimulus_up));repmat(y_min,1,numel(stimulus_down));];
    patch(patch_x,patch_y,'r','FaceAlpha',.2,'EdgeColor','none');  
    ylabel('Displacement (A.U.)');
    if count == numel(rp)
        xlabel('Frames (TR)');
    end
    set(gca,'FontSize',8,'FontWeight','bold');
    set(gca, 'LineWidth',1.5);
    axis([0 numel(stimulus) y_min y_max])
    box off;
    h = legend(f,legend_txt{count});
    pos = 0.20 + 0.12*(count-1);
    set(h,'Position',[pos 0.03 0.070 0.05]);
    set(h,'Box','off');

end
text(numel(Time_series)*10/11,0.07,'(Mean¡ÀSEM)','FontWeight','bold')



end