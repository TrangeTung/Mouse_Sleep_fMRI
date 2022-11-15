
function MY_plot_stimulus_related_signal_variance(rp,stimulus,color)
%% rp
for count = 1:numel(rp)
    TimeCourse_mean = rp{count}; 
    number = size(TimeCourse_mean,2);
    Time_series = mean(TimeCourse_mean,2);
    upper = Time_series+std(TimeCourse_mean,0,2)/sqrt(number);
    lower = flipud(Time_series-std(TimeCourse_mean,0,2)/sqrt(number));
    y_min = 200;
    y_max = 400;%0.6;1.5
    
%     subplot(numel(rp),1,count);
    f = plot([1:numel(Time_series)],Time_series,'color',color/255,'LineWidth',1);
    x = [1:numel(Time_series),fliplr(1:numel(Time_series))]';
    y = [upper;lower];
    patch(x,y,color/255,'FaceAlpha',.1,'EdgeColor','none');  
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
    patch(patch_x,patch_y,'k','FaceAlpha',.2,'EdgeColor','none');  
    set(gca,'FontSize',8,'FontWeight','bold');
    set(gca, 'LineWidth',1.5);
    axis([0 numel(stimulus) y_min y_max])
    box off;
end
end